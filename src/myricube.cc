// Myricube: Experimental voxel renderer using hybrid raycast/mesh
// OpenGL rendering.
//
// David Zhao Akeley 2020 (not a great year)
//
// This was inspired by my discontent with Minecraft's low render
// distance, even on high settings.
//
// Like most voxel renderers, this renderer subdivides an "infinite"
// grid of voxels into cubic chunks. Nearby chunks are drawn using
// conventional methods (by converting the voxels into a mesh of
// triangles). Farther chunks are instead drawn using raycasting: a
// minimal AABB is calculated that contains all the solid voxels of
// the chunk, and the AABB itself is drawn using a raycasting fragment
// shader that checks for collisions with the voxels contained in the
// AABB. The voxel data itself (when raycasting is used) is stored
// using a 3D texture^ -- this is considerably more memory efficient
// than a triangle mesh. (Actually, now that I switched to instanced
// rendering, this memory point may no longer be true...)
//
//     ^ Really, a 3D array in an SSBO; much faster to write than textures.
//
// Cubes of chunks are organized into larger chunk groups. For OpenGL
// efficiency, chunks within the same chunk group share GPU resources
// (3D textures, buffer storage, etc.). Points in space are frequently
// expressed as "group" and "residue" coordinates -- these are the
// coordinates (in chunk-group-size units) of the lower-left of the
// group the point is in, and the remaining offset (in voxel-size
// units) within the group. (See split_coordinate for the precise
// floor-based definition). For example, if chunk groups were 100 x
// 100 x 100 voxels, then point (1, 502.5, -1) has group coordinate
// (0, 5, -1) and residue (1, 2.5, 99).
//
// To be precise, note also that a voxel at coordinate (x,y,z) occupies
// the cube from (x,y,z) to (x+1,y+1,z+1) in space.
//
// The primary benefits of raycasting are the aforementioned memory
// efficiency and the reduced number of triangle vertices processed
// (12 triangles are shared for all voxels per chunk, instead of up to
// 12 triangles per solid voxel). There are several drawbacks that make
// raycasting suitable only for distant chunks:
//
// * Raycasting cost is roughly linear in the size of the rendered
//   triangles (because the raycast fragment shader is the most
//   expensive part). This makes raycasting very expensive for large
//   nearby triangles.
//
// * Raycasting is more prone to rounding errors. This creates
//   unsightly gaps between voxels that are more obvious when nearby.
//
// * Raycasting causes "incorrect" writes to the depth buffer (because
//   the depth by default is that of the AABB, and not the actual
//   voxel), so raycast voxels may incorrectly occlude other geometry
//   (imagine a Minecraft mob standing on a voxel in the center of a
//   chunk). There are potential solutions (e.g. manual gl_FragDepth)
//   but the cost in lost optimizations is considerable.
//
// * Circular reasoning, but in general I've designed the raycaster to
//   be cheap-but-ugly.
//
// This motivates the hybrid renderer I've come up with.

#include "myricube.hh"

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <string.h>
#include <utility>
#include <vector>

#include "app.hh"
#include "camera.hh"
#include "chunk.hh"
#include "hexload.hh"
#include "renderer.hh"
#include "window.hh"

namespace myricube {

// Absolute path of the executable, minus the -bin or .exe, plus -data/
// This is where shaders and stuff are stored.
std::string data_directory;

std::string expand_filename(const std::string& in)
{
    if (data_directory.size() == 0) {
        throw std::logic_error("Cannot call expand_filename before main");
    }
    return in[0] == '/' ? in : data_directory + in;
}

bool ends_with_bin_or_exe(const std::string& in)
{
    auto sz = in.size();
    if (sz < 4) return false;
    const char* suffix = &in[sz - 4];
    return strcmp(suffix, "-bin") == 0 or strcmp(suffix, ".exe") == 0;
}

bool paused = false;

void add_key_targets(Window& window, Camera& camera)
{
    static float speed = 8.0f;
    static float sprint_mod = 1.0f;

    struct Position
    {
        glm::dvec3 eye = glm::dvec3(0);
        float theta = 1.5707f;
        float phi = 1.5707f;
    };
    static Position old_positions_ring_buffer[256];
    static Position future_positions_ring_buffer[256];
    static uint8_t old_idx = 0;
    static uint8_t future_idx = 0;

    static auto get_camera_position = [&camera] () -> Position
    {
        Position p;
        p.eye = camera.get_eye();
        p.theta = camera.get_theta();
        p.phi = camera.get_phi();
        return p;
    };

    static auto push_camera_position = [&]
    {
        old_positions_ring_buffer[--old_idx] = get_camera_position();
    };

    static auto push_camera_position_callback = [&] (KeyArg arg)
    {
        if (arg.repeat) return false;
        push_camera_position();
        return true;
    };

    KeyTarget pop_old_camera, pop_future_camera;
    pop_old_camera.down = [&] (KeyArg) -> bool
    {
        future_positions_ring_buffer[--future_idx] =
            get_camera_position();
        Position p = old_positions_ring_buffer[old_idx++];
        camera.set_eye(p.eye);
        camera.set_theta(p.theta);
        camera.set_phi(p.phi);
        return true;
    };
    pop_future_camera.down = [&] (KeyArg) -> bool
    {
        old_positions_ring_buffer[--old_idx] =
            get_camera_position();
        Position p = future_positions_ring_buffer[future_idx++];
        camera.set_eye(p.eye);
        camera.set_theta(p.theta);
        camera.set_phi(p.phi);
        return true;
    };
    window.add_key_target("pop_old_camera", pop_old_camera);
    window.add_key_target("pop_future_camera", pop_future_camera);

    KeyTarget forward, backward, leftward, rightward, upward, downward;
    forward.down = push_camera_position_callback;
    forward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera.frenet_move(0, 0, +arg.dt * speed * sprint_mod);
        return true;
    };
    backward.down = push_camera_position_callback;
    backward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera.frenet_move(0, 0, -arg.dt * speed * sprint_mod);
        return true;
    };
    leftward.down = push_camera_position_callback;
    leftward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera.frenet_move(-arg.dt * speed * sprint_mod, 0, 0);
        return true;
    };
    rightward.down = push_camera_position_callback;
    rightward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera.frenet_move(+arg.dt * speed * sprint_mod, 0, 0);
        return true;
    };
    upward.down = push_camera_position_callback;
    upward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera.frenet_move(0, +arg.dt * speed * sprint_mod, 0);
        return true;
    };
    downward.down = push_camera_position_callback;
    downward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera.frenet_move(0, -arg.dt * speed * sprint_mod, 0);
        return true;
    };
    window.add_key_target("forward", forward);
    window.add_key_target("backward", backward);
    window.add_key_target("leftward", leftward);
    window.add_key_target("rightward", rightward);
    window.add_key_target("upward", upward);
    window.add_key_target("downward", downward);

    KeyTarget sprint, speed_up, slow_down;
    sprint.down = [&] (KeyArg) -> bool
    {
        sprint_mod = 7.0f;
        return true;
    };
    sprint.up = [&] (KeyArg) -> bool
    {
        sprint_mod = 1.0f;
        return true;
    };
    speed_up.down = [&] (KeyArg arg) -> bool
    {
        if (!arg.repeat) speed *= 2.0f;
        return !arg.repeat;
    };
    slow_down.down = [&] (KeyArg arg) -> bool
    {
        if (!arg.repeat) speed *= 0.5f;
        return !arg.repeat;
    };
    window.add_key_target("sprint", sprint);
    window.add_key_target("speed_up", speed_up);
    window.add_key_target("slow_down", slow_down);

    KeyTarget vertical_scroll, horizontal_scroll, look_around;
    look_around.down = push_camera_position_callback;
    look_around.per_frame = [&] (KeyArg arg) -> bool
    {
        camera.inc_theta(arg.mouse_rel_x * arg.dt * 0.01f);
        camera.inc_phi(arg.mouse_rel_y * arg.dt * 0.01f);
        return true;
    };
    vertical_scroll.down = [&] (KeyArg arg) -> bool
    {
        camera.inc_phi(arg.amount * -0.05f);
        return true;
    };
    horizontal_scroll.down = [&] (KeyArg arg) -> bool
    {
        camera.inc_theta(arg.amount * -0.05f);
        return true;
    };
    window.add_key_target("look_around", look_around);
    window.add_key_target("vertical_scroll", vertical_scroll);
    window.add_key_target("horizontal_scroll", horizontal_scroll);

    KeyTarget pause;
    pause.down = [&] (KeyArg) -> bool
    {
        paused = !paused;
        return true;
    };
    window.add_key_target("pause", pause);

    KeyTarget toggle_fog;
    toggle_fog.down = [&camera] (KeyArg) -> bool
    {
        camera.set_fog(!camera.get_fog());
        return true;
    };
    window.add_key_target("toggle_fog", toggle_fog);

    KeyTarget toggle_black_fog;
    toggle_black_fog.down = [&camera] (KeyArg) -> bool
    {
        camera.use_black_fog(!camera.use_black_fog());
        return true;
    };
    window.add_key_target("toggle_black_fog", toggle_black_fog);

    // Maybe I should dehackify this variable one day.
    extern bool chunk_debug;
    KeyTarget toggle_chunk_debug;
    toggle_chunk_debug.down = [&] (KeyArg) -> bool
    {
        chunk_debug = !chunk_debug;
        return true;
    };
    window.add_key_target("toggle_chunk_debug", toggle_chunk_debug);

    extern bool evict_stats_debug;
    KeyTarget toggle_evict_stats_debug;
    toggle_evict_stats_debug.down = [&] (KeyArg) -> bool
    {
        evict_stats_debug = !evict_stats_debug;
        return true;
    };
    window.add_key_target("toggle_evict_stats_debug", toggle_evict_stats_debug);

    KeyTarget unload;
    unload.down = [&] (KeyArg) -> bool
    {
        camera.unload_gpu_storage();
        evict_stats_debug = true;
        return true;
    };
    window.add_key_target("unload_gpu_storage", unload);

    KeyTarget increase_far_plane, decrease_far_plane;
    increase_far_plane.down = [&] (KeyArg arg) -> bool
    {
        auto new_far_plane = 64 + camera.get_far_plane();

        if (!arg.repeat) {
            camera.set_far_plane(new_far_plane);
            fprintf(stderr, "Far plane: %i\n", int(new_far_plane));
        }
        return !arg.repeat;
    };
    window.add_key_target("increase_far_plane", increase_far_plane);

    decrease_far_plane.down = [&] (KeyArg arg) -> bool
    {
        auto new_far_plane = -64 + camera.get_far_plane();
        if (new_far_plane < 64) new_far_plane = 64;

        if (!arg.repeat) {
            camera.set_far_plane(new_far_plane);
            fprintf(stderr, "Far plane: %i\n", int(new_far_plane));
        }
        return !arg.repeat;
    };
    window.add_key_target("decrease_far_plane", decrease_far_plane);
}

// Given the full path of a key binds file, parse it for key bindings
// and add it to the window's database of key bindings (physical
// key/mouse button to KeyTarget name associations).
//
// Syntax: the file should consist of lines of pairs of key names and
// KeyTarget names. Blank (all whitespace) lines are allowed as well
// as comments, which go from a # character to the end of the line.
//
// Returns true iff successful (check errno on false).
bool add_key_binds_from_file(Window& window, std::string filename) noexcept
{
    FILE* file = fopen(filename.c_str(), "r");
    if (file == nullptr) {
        fprintf(stderr, "Could not open %s\n", filename.c_str());
        return false;
    }

    int line_number = 0;

    auto skip_whitespace = [file]
    {
        int c;
        while (1) {
            c = fgetc(file);
            if (c == EOF) return;
            if (c == '\n' or !isspace(c)) {
                ungetc(c, file);
                return;
            }
        }
    };
    errno = 0;

    bool eof = false;
    while (!eof) {
        std::string key_name;
        std::string target_name;
        ++line_number;

        int c;
        skip_whitespace();

        // Parse key name (not case sensitive -- converted to lower case)
        while (1) {
            c = fgetc(file);

            if (c == EOF) {
                if (errno != 0 and errno != EAGAIN) goto bad_eof;
                eof = true;
                goto end_line;
            }
            if (c == '\n') goto end_line;
            if (isspace(c)) break;
            if (c == '#') goto comment;
            key_name.push_back(c);
        }

        skip_whitespace();

        // Parse target name (case sensitive)
        while (1) {
            c = fgetc(file);

            if (c == EOF) {
                if (errno != 0 and errno != EAGAIN) goto bad_eof;
                eof = true;
                goto end_line;
            }
            if (c == '\n') goto end_line;
            if (isspace(c)) break;
            if (c == '#') goto comment;
            target_name.push_back(c);
        }

        skip_whitespace();

        // Check for unexpected cruft at end of line.
        c = fgetc(file);
        if (c == EOF) {
            if (errno != 0 and errno != EAGAIN) goto bad_eof;
            eof = true;
            goto end_line;
        }
        else if (c == '#') {
            goto comment;
        }
        else if (c == '\n') {
            goto end_line;
        }
        else {
            fprintf(stderr, "%s:%i unexpected third token"
                " starting with '%c'\n",
                filename.c_str(), line_number, c);
            errno = EINVAL;
            goto bad_eof;
        }

        // Skip over comment characters from # to \n
      comment:
        while (1) {
            c = fgetc(file);
            if (c == EOF) {
                if (errno != 0 and errno != EAGAIN) goto bad_eof;
                eof = true;
                goto end_line;
            }
            if (c == '\n') {
                break;
            }
        }
      end_line:
        // skip blank lines silently.
        if (key_name.size() == 0) continue;

        // Complain if only one token is provided on a line.
        if (target_name.size() == 0) {
            fprintf(stderr, "%s:%i key name without target name.\n",
                filename.c_str(), line_number);
            errno = EINVAL;
            goto bad_eof;
        }

        auto keycode = keycode_from_name(key_name);
        if (keycode == 0) {
            fprintf(stderr, "%s:%i unknown key name %s.\n",
                filename.c_str(), line_number, key_name.c_str());
            errno = EINVAL;
            goto bad_eof;
        }

        fprintf(stderr, "Binding %s (%i) to %s\n",
            key_name.c_str(), keycode, target_name.c_str());
        window.bind_keycode(keycode, target_name);
    }

    if (fclose(file) != 0) {
        fprintf(stderr, "Error closing %s\n", filename.c_str());
        return false;
    }
    return true;
  bad_eof:
    fprintf(stderr, "Warning: unexpected end of parsing.\n");
    int eof_errno = errno;
    fclose(file);
    errno = eof_errno;
    return true; // I'm getting bogus EOF fails all the time so fake success :/
}

void bind_keys(Window& window)
{
    auto default_file = expand_filename("default-keybinds.txt");
    auto user_file = expand_filename("keybinds.txt");

    bool default_okay = add_key_binds_from_file(window, default_file);
    if (!default_okay) {
        fprintf(stderr, "Failed to parse %s\n", default_file.c_str());
        fprintf(stderr, "%s (%i)\n", strerror(errno), errno);
        exit(2);
    }

    bool user_okay = add_key_binds_from_file(window, user_file);
    if (!user_okay) {
        if (errno == ENOENT) {
            fprintf(stderr, "Custom keybinds file %s not found.\n",
                user_file.c_str());
        }
        else {
            fprintf(stderr, "Failed to parse %s\n", user_file.c_str());
            fprintf(stderr, "%s (%i)\n", strerror(errno), errno);
            exit(2);
        }
    }
}

void set_window_title(Window& window)
{
    const char format[] = "Myricube %07.2f FPS %03dms frame time";
    char title[sizeof format + 20];
    sprintf(title, format, window.get_fps(), window.get_frame_time_ms());
    window.set_title(title);
}

int Main(std::vector<std::string> args)
{
    if (args.at(0)[0] != '/') {
        fprintf(stderr, "Warning: %s expected to be absolute path\n"
            "(call through wrapper script).\n", args[0].c_str());
    }
    // Data directory (where shaders are stored) is the path of this
    // executable, with the -bin or .exe file extension replaced with
    // -data. Construct that directory name here.
    data_directory = args[0];
    if (!ends_with_bin_or_exe(data_directory)) {
        fprintf(stderr, "%s should end with '-bin' or '.exe'\n",
            args[0].c_str());
        return 1;
    }
    for (int i = 0; i < 4; ++i) data_directory.pop_back();
    data_directory += "-data/";

    // Instantiate the camera.
    Camera camera;

    // Create a window; callback ensures these window dimensions stay accurate.
    int screen_x = 0, screen_y = 0;
    auto on_window_resize = [&camera, &screen_x, &screen_y] (int x, int y)
    {
        viewport(x, y);
        camera.set_window_size(x, y);
        screen_x = x;
        screen_y = y;
    };
    Window window(on_window_resize);

    // Set up keyboard controls.
    add_key_targets(window, camera);
    bind_keys(window);

    // Instantiate the app (which owns the VoxelWorld to render) based
    // on the user's myricube_app environment variable, if any.
    const char* app_name = getenv("myricube_app");
    if (app_name == nullptr) app_name = "RandomWalk";
    std::unique_ptr<App> app(new_named_app(app_name));
    if (app == nullptr) {
        fprintf(stderr, "Known app names:\n");
        stderr_dump_app_names();
        panic("myricube_app environment variable set to unknown name",
            app_name);
    }
    app->add_key_targets(window);
    float dt = 0;
    VoxelWorld* world = &app->update(dt);

    // Render loop.
    gl_first_time_setup();
    while (window.update_swap_buffers(&dt)) {
        if (!paused) world = &app->update(dt);
        camera.fix_dirty();

        gl_clear();
        bind_global_f32_depth_framebuffer(screen_x, screen_y);
        gl_clear();
        render_world_mesh_step(*world, camera);
        render_world_raycast_step(*world, camera);
        render_background(camera);
        finish_global_f32_depth_framebuffer(screen_x, screen_y);

        extern bool evict_stats_debug;
        evict_stats_debug = false;

        set_window_title(window);
    }

    // Test out the silly hexload class.
    window.set_title("Autosaving...");
    auto autosave_filename = expand_filename("autosave.myricube.hex");
    bool okay = write_hex(*world, autosave_filename);
    if (!okay) {
        fprintf(stderr, "Too bad, autosave failed.\n");
    }

    return 0;
}

} // end namespace

int main(int argc, char** argv)
{
    std::vector<std::string> args;
    for (int i = 0; i < argc; ++i) {
        args.emplace_back(argv[i]);
    }
    return myricube::Main(std::move(args));
}
