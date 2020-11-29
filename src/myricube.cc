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
// using a 3D texture.
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
// The primary benefit of raycasting is the reduced number of triangle
// vertices processed (12 triangles are shared for all voxels per
// chunk, instead of up to 12 triangles per solid voxel). There are
// several drawbacks that make raycasting suitable only for distant
// chunks:
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
#include "RendererGL.hh"
#include "RendererVk.hh"
#include "RenderThread.hh"
#include "window.hh"
#include "voxels.hh"

namespace myricube {

// Absolute path of the executable, minus the -bin or .exe, plus -data/
// This is where shaders and stuff are stored.
filename_string data_directory;

filename_string expand_filename(const std::string& in)
{
    if (data_directory.size() == 0) {
        throw std::logic_error("Cannot call expand_filename before main");
    }
    if (in.empty()) {
        throw std::runtime_error("Empty expand_filename argument");
    }
    if (in[0] == '/') {
        throw std::runtime_error("Forbid absolute path in expand_filename");
    }
    return filename_concat_c_str(data_directory, in.c_str());
}

bool ends_with_bin_or_exe(const filename_string& in)
{
    auto sz = in.size();
    if (sz < 4) return false;
    const filename_char* suffix = &in[sz - 4];
#ifdef MYRICUBE_WINDOWS
    return wcscmp(suffix, L"-bin") == 0 or wcscmp(suffix, L".exe") == 0;
#else
    return strcmp(suffix, "-bin") == 0 or strcmp(suffix, ".exe") == 0;
#endif
}

bool paused = false;

void add_key_targets(Window& window, std::shared_ptr<SyncCamera> camera_arg)
{
    static std::shared_ptr<SyncCamera> camera = std::move(camera_arg);
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

    static auto get_camera_position = [&] () -> Position
    {
        Position p;
        p.eye = camera->get_eye();
        p.theta = camera->get_theta();
        p.phi = camera->get_phi();
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
        camera->set_eye(p.eye);
        camera->set_theta(p.theta);
        camera->set_phi(p.phi);
        return true;
    };
    pop_future_camera.down = [&] (KeyArg) -> bool
    {
        old_positions_ring_buffer[--old_idx] =
            get_camera_position();
        Position p = future_positions_ring_buffer[future_idx++];
        camera->set_eye(p.eye);
        camera->set_theta(p.theta);
        camera->set_phi(p.phi);
        return true;
    };
    window.add_key_target("pop_old_camera", pop_old_camera);
    window.add_key_target("pop_future_camera", pop_future_camera);

    KeyTarget forward, backward, leftward, rightward, upward, downward;
    forward.down = push_camera_position_callback;
    forward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera->frenet_move(0, 0, +arg.dt * speed * sprint_mod);
        return true;
    };
    backward.down = push_camera_position_callback;
    backward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera->frenet_move(0, 0, -arg.dt * speed * sprint_mod);
        return true;
    };
    leftward.down = push_camera_position_callback;
    leftward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera->frenet_move(-arg.dt * speed * sprint_mod, 0, 0);
        return true;
    };
    rightward.down = push_camera_position_callback;
    rightward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera->frenet_move(+arg.dt * speed * sprint_mod, 0, 0);
        return true;
    };
    upward.down = push_camera_position_callback;
    upward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera->frenet_move(0, +arg.dt * speed * sprint_mod, 0);
        return true;
    };
    downward.down = push_camera_position_callback;
    downward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera->frenet_move(0, -arg.dt * speed * sprint_mod, 0);
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
        camera->inc_theta(arg.mouse_rel_x * arg.dt * 0.01f);
        camera->inc_phi(arg.mouse_rel_y * arg.dt * 0.01f);
        return true;
    };
    vertical_scroll.down = [&] (KeyArg arg) -> bool
    {
        camera->inc_phi(arg.amount * -0.05f);
        return true;
    };
    horizontal_scroll.down = [&] (KeyArg arg) -> bool
    {
        camera->inc_theta(arg.amount * -0.05f);
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
    toggle_fog.down = [&] (KeyArg) -> bool
    {
        camera->set_fog(!camera->get_fog());
        return true;
    };
    window.add_key_target("toggle_fog", toggle_fog);

    KeyTarget toggle_black_fog;
    toggle_black_fog.down = [&] (KeyArg) -> bool
    {
        camera->use_black_fog(!camera->use_black_fog());
        return true;
    };
    window.add_key_target("toggle_black_fog", toggle_black_fog);

    KeyTarget toggle_chunk_debug;
    toggle_chunk_debug.down = [&] (KeyArg) -> bool
    {
        camera->set_chunk_debug(!camera->is_chunk_debug());
        return true;
    };
    window.add_key_target("toggle_chunk_debug", toggle_chunk_debug);

    KeyTarget increase_far_plane, decrease_far_plane;
    increase_far_plane.down = [&] (KeyArg arg) -> bool
    {
        auto new_far_plane = 64 + camera->get_far_plane();

        if (!arg.repeat) {
            camera->set_far_plane(new_far_plane);
            fprintf(stderr, "Far plane: %i\n", int(new_far_plane));
        }
        return !arg.repeat;
    };
    window.add_key_target("increase_far_plane", increase_far_plane);

    decrease_far_plane.down = [&] (KeyArg arg) -> bool
    {
        auto new_far_plane = -64 + camera->get_far_plane();
        if (new_far_plane < 64) new_far_plane = 64;

        if (!arg.repeat) {
            camera->set_far_plane(new_far_plane);
            fprintf(stderr, "Far plane: %i\n", int(new_far_plane));
        }
        return !arg.repeat;
    };
    window.add_key_target("decrease_far_plane", decrease_far_plane);

    KeyTarget increase_target_fragments, decrease_target_fragments;
    constexpr int min_nonzero_fragments = 250000;

    increase_target_fragments.down = [&] (KeyArg) -> bool
    {
        if (camera->get_target_fragments() <= 0) {
            camera->set_target_fragments(min_nonzero_fragments);
        }
        else camera->set_target_fragments(camera->get_target_fragments() * 2);

        fprintf(stderr, "%i target fragments.\n", camera->get_target_fragments());
        return true;
    };
    window.add_key_target(
        "increase_target_fragments", increase_target_fragments);
    decrease_target_fragments.down = [&] (KeyArg) -> bool
    {
        if (camera->get_target_fragments() <= min_nonzero_fragments) {
            camera->set_target_fragments(0);
            fprintf(stderr, "Unlimited fragments.\n");
        }
        else {
            camera->set_target_fragments(camera->get_target_fragments() / 2);
            fprintf(stderr, "%i target fragments.\n", camera->get_target_fragments());
        }
        return true;
    };
    window.add_key_target(
        "decrease_target_fragments", decrease_target_fragments);
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
bool add_key_binds_from_file(Window& window, filename_string filename) noexcept
{
// Hack for printing 8/16-bit strings on Windows, 8-bit only on Linux.
// Open with 16-bit filename on Windows, 8-bit on Linux.
#ifdef MYRICUBE_WINDOWS
    #define WS "%ls"
    #define HS "%s"
    FILE* file = _wfopen(filename.c_str(), L"r");
#else
    #define WS "%s"
    #define HS "%s"
    FILE* file = fopen(filename.c_str(), "r");
#endif

    if (file == nullptr) {
        fprintf(stderr, "Could not open " WS "\n", filename.c_str());
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
            fprintf(stderr, WS ":%i unexpected third token"
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
            fprintf(stderr, WS ":%i key name without target name.\n",
                filename.c_str(), line_number);
            errno = EINVAL;
            goto bad_eof;
        }

        auto keycode = keycode_from_name(key_name);
        if (keycode == 0) {
            fprintf(stderr, WS ":%i unknown key name " HS ".\n",
                filename.c_str(), line_number, key_name.c_str());
            errno = EINVAL;
            goto bad_eof;
        }

        fprintf(stderr, "Binding " HS " (%i) to " HS "\n",
            key_name.c_str(), keycode, target_name.c_str());
        window.bind_keycode(keycode, target_name);
    }

    if (fclose(file) != 0) {
        fprintf(stderr, "Error closing " WS "\n", filename.c_str());
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
    filename_string default_file = expand_filename("default-keybinds.txt");
    filename_string user_file = expand_filename("keybinds.txt");

    bool default_okay = add_key_binds_from_file(window, default_file);
    if (!default_okay) {
        fprintf(stderr, "Failed to parse " WS "\n", default_file.c_str());
        fprintf(stderr, "%s (%i)\n", strerror(errno), errno);
        exit(2);
    }

    bool user_okay = add_key_binds_from_file(window, user_file);
    if (!user_okay) {
        if (errno == ENOENT) {
            fprintf(stderr, "Custom keybinds file " WS " not found.\n",
                user_file.c_str());
        }
        else {
            fprintf(stderr, "Failed to parse " WS "\n", user_file.c_str());
            fprintf(stderr, "%s (%i)\n", strerror(errno), errno);
            exit(2);
        }
    }
}

void set_window_title(Window& window, const RenderThread& renderer)
{
    const char format[] = "Myricube %07.2f FPS %03dms frame time";
    char title[sizeof format + 20];
    int milliseconds = int(renderer.get_frame_time() * 1000);
    sprintf(title, format, renderer.get_fps(), milliseconds);
    window.set_title(title);
}

int Main(std::vector<filename_string> args)
{
#ifndef MYRICUBE_WINDOWS
    if (args.at(0)[0] != '/') {
        fprintf(stderr, "Warning: %s expected to be absolute path\n"
            "(call through wrapper script).\n", args[0].c_str());
    }
#endif
    // Data directory (where shaders are stored) is the path of this
    // executable, with the -bin or .exe file extension replaced with
    // -data. Construct that directory name here.
    data_directory = args[0];
    if (!ends_with_bin_or_exe(data_directory)) {
        fprintf(stderr, WS " should end with '-bin' or '.exe'\n",
            args[0].c_str());
        return 1;
    }
    for (int i = 0; i < 4; ++i) data_directory.pop_back();
    data_directory = filename_concat_c_str(data_directory, "-data/");

    // Check if we are using OpenGL.
    bool use_OpenGL = false;
    const char* api_name = getenv("myricube_api");
    if (api_name != nullptr) {
        if (strcmp(api_name, "gl") == 0) {
            use_OpenGL = true;
        }
        else if (strcmp(api_name, "vk") != 0) {
            panic("myricube_api environment variable set to unknown value",
                api_name);
        }
    }

    // Instantiate the shared camera (shared with the renderer).
    std::shared_ptr<SyncCamera> camera_ptr(new SyncCamera());

    // Create a window; callback ensures the window dimensions stay accurate.
    auto on_window_resize = [&camera_ptr] (int x, int y)
    {
        camera_ptr->set_framebuffer_size(x, y);
    };
    std::shared_ptr<Window> window_ptr(
        new Window(on_window_resize, use_OpenGL));
    Window& window = *window_ptr;

    // Set up keyboard controls.
    add_key_targets(window, camera_ptr);
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
    WorldHandle handle = world->get_handle();
    filename_string expected_world_dir = handle.get_directory_name();

    // Start renderer thread.
    RenderArgs render_args { window_ptr, camera_ptr, handle };
    RenderThread renderer(
        use_OpenGL ? RendererGL_Factory : RendererVk_Factory, render_args);
    camera_ptr->set_max_frame_new_chunk_groups(32);

    // Meanwhile main thread handles user input and runs the demonstration app.
    while (window.update_events(&dt)) {
        if (!paused) world = &app->update(dt);

        if (world->get_handle().get_directory_name() != expected_world_dir) {
            throw std::runtime_error("Not implemented: "
                "App changing voxel world directory");
        }

        set_window_title(window, renderer);
    }

    // Stop renderer thread (RenderThread destructor implicit).
    return 0;
}

} // end namespace

#ifdef MYRICUBE_WINDOWS
int WinMain(int argc, wchar_t** argv)
{
    std::vector<myricube::filename_string> args;
    for (int i = 0; i < argc; ++i) {
        args.emplace_back(argv[i]);
    }
    return myricube::Main(std::move(args));
}
#else
int main(int argc, char** argv)
{
    std::vector<std::string> args;
    for (int i = 0; i < argc; ++i) {
        args.emplace_back(argv[i]);
    }
    return myricube::Main(std::move(args));
}
#endif
