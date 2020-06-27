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
// using a 3D texture -- this is considerably more memory efficient
// than a triangle mesh.
//
// Cubes of chunks are organized into larger chunk groups. For OpenGL
// efficiency, chunks within the same chunk group share GPU resources
// (3D textures, buffer storage, etc.). Points in space are frequently
// expressed as "group" and "residue" coordinates -- these are the
// coordinates (in chunk-group-size units) of the lower-left of the
// chunk the point is in, and the remaining offset within the
// chunk. (See split_coordinate for the precise floor-based
// definition). For example, if chunk groups were 100 x 100 x 100
// voxels, then point (1, 502.5, -1) has chunk coordinate (0, 5, -1)
// and residue (1, 2.5, 99).
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

#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "app.hh"
#include "camera.hh"
#include "chunk.hh"
#include "renderer.hh"
#include "window.hh"

#include "SDL2/SDL.h"

namespace myricube {

// Absolute path of the executable, minus the -bin, plus -data/
// This is where shaders and stuff are stored.
std::string data_directory;

std::string expand_filename(const std::string& in)
{
    if (data_directory.size() == 0) {
        throw std::logic_error("Cannot call expand_filename before main");
    }
    return in[0] == '/' ? in : data_directory + in;
}

bool ends_with_dash_bin(const std::string& in)
{
    auto sz = in.size();
    return sz >= 4 and
           in[sz-4] == '-' and
           in[sz-3] == 'b' and
           in[sz-2] == 'i' and
           in[sz-1] == 'n';
           // sigh...
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

    extern bool chunk_debug;
    KeyTarget toggle_chunk_debug;
    toggle_chunk_debug.down = [&] (KeyArg) -> bool
    {
        chunk_debug = !chunk_debug;
        return true;
    };
    window.add_key_target("toggle_chunk_debug", toggle_chunk_debug);

    KeyTarget toggle_culling_freeze_target;
    toggle_culling_freeze_target.down = [&] (KeyArg) -> bool
    {
        toggle_culling_freeze(camera);
        return true;
    };
    window.add_key_target("toggle_culling_freeze", toggle_culling_freeze_target);

    KeyTarget unload;
    unload.down = [&] (KeyArg) -> bool
    {
        camera.unload_gpu_storage();
        return true;
    };
    window.add_key_target("unload_gpu_storage", unload);
}

extern std::unordered_map<std::string, int> key_name_to_key_code_map;

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
            key_name.push_back(tolower(c));
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

        auto it = key_name_to_key_code_map.find(key_name);
        if (it == key_name_to_key_code_map.end()) {
            fprintf(stderr, "%s:%i unknown key name %s.\n",
                filename.c_str(), line_number, key_name.c_str());
            errno = EINVAL;
            goto bad_eof;
        }

        fprintf(stderr, "Binding %s (%i) to %s\n",
            key_name.c_str(), it->second, target_name.c_str());
        window.bind_keycode(it->second, target_name);
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
    std::string title = "Myricube ";
    title += std::to_string(window.get_fps()) + " FPS ";
    title += std::to_string(window.get_frame_time_ms()) + "ms frame time";

    window.set_title(title);
}

int Main(std::vector<std::string> args)
{
    if (args.at(0)[0] != '/') {
        fprintf(stderr, "%s should be absolute path\n"
            "(call through wrapper script).\n", args[0].c_str());
        return 1;
    }
    data_directory = args[0];
    // if (!data_directory.ends_with("-bin")) {
    if (!ends_with_dash_bin(data_directory)) {
        fprintf(stderr, "%s should end with '-bin'\n",
            args[0].c_str());
        return 1;
    }
    for (int i = 0; i < 4; ++i) data_directory.pop_back();
    data_directory += "-data/";

    VoxelWorld world;
    Camera camera;

    auto on_window_resize = [&camera] (int x, int y)
    {
        viewport(x, y);
        camera.set_window_size(x, y);
    };
    Window window(on_window_resize);
    add_key_targets(window, camera);
    bind_keys(window);

    app_init(world, window);
    gl_first_time_setup();
    while (window.update_swap_buffers(5)) {
        if (!paused) app_update(world);
        gl_clear();
        camera.fix_dirty();
        render_world_mesh_step(world, camera);
        render_world_raycast_step(world, camera);
        set_window_title(window);
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

namespace myricube {

std::unordered_map<std::string, int> key_name_to_key_code_map {
    { "mouse-1", -1 },
    { "left-mouse", -1 },
    { "mouse-2", -2 },
    { "middle-mouse", -2 },
    { "mouse-3", -3 },
    { "right-mouse", -3 },
    { "mouse-4", -4 },
    { "scroll-up", -4 },
    { "mouse-5", -5 },
    { "scroll-down", -5 },
    { "mouse-6", -6 },
    { "scroll-left", -6 },
    { "mouse-7", -7 },
    { "scroll-right", -7 },
    { "mouse-8", -8 },
    { "thumb-button", -8 },
    { "x1", -8 },
    { "mouse-9", -9 },
    { "thumb-button-2", -9 },
    { "x2", -9 },
    { "mouse-10", -10 },

    { "a", SDL_SCANCODE_A },
    { "b", SDL_SCANCODE_B },
    { "c", SDL_SCANCODE_C },
    { "d", SDL_SCANCODE_D },
    { "e", SDL_SCANCODE_E },
    { "f", SDL_SCANCODE_F },
    { "g", SDL_SCANCODE_G },
    { "h", SDL_SCANCODE_H },
    { "i", SDL_SCANCODE_I },
    { "j", SDL_SCANCODE_J },
    { "k", SDL_SCANCODE_K },
    { "l", SDL_SCANCODE_L },
    { "m", SDL_SCANCODE_M },
    { "n", SDL_SCANCODE_N },
    { "o", SDL_SCANCODE_O },
    { "p", SDL_SCANCODE_P },
    { "q", SDL_SCANCODE_Q },
    { "r", SDL_SCANCODE_R },
    { "s", SDL_SCANCODE_S },
    { "t", SDL_SCANCODE_T },
    { "u", SDL_SCANCODE_U },
    { "v", SDL_SCANCODE_V },
    { "w", SDL_SCANCODE_W },
    { "x", SDL_SCANCODE_X },
    { "y", SDL_SCANCODE_Y },
    { "z", SDL_SCANCODE_Z },

    { "1", SDL_SCANCODE_1 },
    { "2", SDL_SCANCODE_2 },
    { "3", SDL_SCANCODE_3 },
    { "4", SDL_SCANCODE_4 },
    { "5", SDL_SCANCODE_5 },
    { "6", SDL_SCANCODE_6 },
    { "7", SDL_SCANCODE_7 },
    { "8", SDL_SCANCODE_8 },
    { "9", SDL_SCANCODE_9 },
    { "0", SDL_SCANCODE_0 },

    { "return", SDL_SCANCODE_RETURN },
    { "escape", SDL_SCANCODE_ESCAPE },
    { "backspace", SDL_SCANCODE_BACKSPACE },
    { "tab", SDL_SCANCODE_TAB },
    { "space", SDL_SCANCODE_SPACE },

    { "minus", SDL_SCANCODE_MINUS },
    { "equals", SDL_SCANCODE_EQUALS },
    { "leftbracket", SDL_SCANCODE_LEFTBRACKET },
    { "rightbracket", SDL_SCANCODE_RIGHTBRACKET },
    { "backslash", SDL_SCANCODE_BACKSLASH },
    { "nonushash", SDL_SCANCODE_NONUSHASH },
    { "semicolon", SDL_SCANCODE_SEMICOLON },
    { "apostrophe", SDL_SCANCODE_APOSTROPHE },
    { "grave", SDL_SCANCODE_GRAVE },
    { "comma", SDL_SCANCODE_COMMA },
    { "period", SDL_SCANCODE_PERIOD },
    { "slash", SDL_SCANCODE_SLASH },

    { "capslock", SDL_SCANCODE_CAPSLOCK },

    { "f1", SDL_SCANCODE_F1 },
    { "f2", SDL_SCANCODE_F2 },
    { "f3", SDL_SCANCODE_F3 },
    { "f4", SDL_SCANCODE_F4 },
    { "f5", SDL_SCANCODE_F5 },
    { "f6", SDL_SCANCODE_F6 },
    { "f7", SDL_SCANCODE_F7 },
    { "f8", SDL_SCANCODE_F8 },
    { "f9", SDL_SCANCODE_F9 },
    { "f10", SDL_SCANCODE_F10 },
    { "f11", SDL_SCANCODE_F11 },
    { "f12", SDL_SCANCODE_F12 },

    { "printscreen", SDL_SCANCODE_PRINTSCREEN },
    { "scrolllock", SDL_SCANCODE_SCROLLLOCK },
    { "pause", SDL_SCANCODE_PAUSE },
    { "insert", SDL_SCANCODE_INSERT },
    { "home", SDL_SCANCODE_HOME },
    { "pageup", SDL_SCANCODE_PAGEUP },
    { "delete", SDL_SCANCODE_DELETE },
    { "end", SDL_SCANCODE_END },
    { "pagedown", SDL_SCANCODE_PAGEDOWN },
    { "right", SDL_SCANCODE_RIGHT },
    { "left", SDL_SCANCODE_LEFT },
    { "down", SDL_SCANCODE_DOWN },
    { "up", SDL_SCANCODE_UP },

    { "numlockclear", SDL_SCANCODE_NUMLOCKCLEAR },
    { "kp_divide", SDL_SCANCODE_KP_DIVIDE },
    { "kp_multiply", SDL_SCANCODE_KP_MULTIPLY },
    { "kp_minus", SDL_SCANCODE_KP_MINUS },
    { "kp_plus", SDL_SCANCODE_KP_PLUS },
    { "kp_enter", SDL_SCANCODE_KP_ENTER },
    { "kp_1", SDL_SCANCODE_KP_1 },
    { "kp_2", SDL_SCANCODE_KP_2 },
    { "kp_3", SDL_SCANCODE_KP_3 },
    { "kp_4", SDL_SCANCODE_KP_4 },
    { "kp_5", SDL_SCANCODE_KP_5 },
    { "kp_6", SDL_SCANCODE_KP_6 },
    { "kp_7", SDL_SCANCODE_KP_7 },
    { "kp_8", SDL_SCANCODE_KP_8 },
    { "kp_9", SDL_SCANCODE_KP_9 },
    { "kp_0", SDL_SCANCODE_KP_0 },
    { "kp_period", SDL_SCANCODE_KP_PERIOD },

    { "nonusbackslash", SDL_SCANCODE_NONUSBACKSLASH },
    { "application", SDL_SCANCODE_APPLICATION },
    { "power", SDL_SCANCODE_POWER },
    { "kp_equals", SDL_SCANCODE_KP_EQUALS },
    { "f13", SDL_SCANCODE_F13 },
    { "f14", SDL_SCANCODE_F14 },
    { "f15", SDL_SCANCODE_F15 },
    { "f16", SDL_SCANCODE_F16 },
    { "f17", SDL_SCANCODE_F17 },
    { "f18", SDL_SCANCODE_F18 },
    { "f19", SDL_SCANCODE_F19 },
    { "f20", SDL_SCANCODE_F20 },
    { "f21", SDL_SCANCODE_F21 },
    { "f22", SDL_SCANCODE_F22 },
    { "f23", SDL_SCANCODE_F23 },
    { "f24", SDL_SCANCODE_F24 },
    { "execute", SDL_SCANCODE_EXECUTE },
    { "help", SDL_SCANCODE_HELP },
    { "menu", SDL_SCANCODE_MENU },
    { "select", SDL_SCANCODE_SELECT },
    { "stop", SDL_SCANCODE_STOP },
    { "again", SDL_SCANCODE_AGAIN },
    { "undo", SDL_SCANCODE_UNDO },
    { "cut", SDL_SCANCODE_CUT },
    { "copy", SDL_SCANCODE_COPY },
    { "paste", SDL_SCANCODE_PASTE },
    { "find", SDL_SCANCODE_FIND },
    { "mute", SDL_SCANCODE_MUTE },
    { "volumeup", SDL_SCANCODE_VOLUMEUP },
    { "volumedown", SDL_SCANCODE_VOLUMEDOWN },
    { "kp_comma", SDL_SCANCODE_KP_COMMA },
    { "kp_equalsas400", SDL_SCANCODE_KP_EQUALSAS400 },

    { "international1", SDL_SCANCODE_INTERNATIONAL1 },
    { "international2", SDL_SCANCODE_INTERNATIONAL2 },
    { "international3", SDL_SCANCODE_INTERNATIONAL3 },
    { "international4", SDL_SCANCODE_INTERNATIONAL4 },
    { "international5", SDL_SCANCODE_INTERNATIONAL5 },
    { "international6", SDL_SCANCODE_INTERNATIONAL6 },
    { "international7", SDL_SCANCODE_INTERNATIONAL7 },
    { "international8", SDL_SCANCODE_INTERNATIONAL8 },
    { "international9", SDL_SCANCODE_INTERNATIONAL9 },
    { "lang1", SDL_SCANCODE_LANG1 },
    { "lang2", SDL_SCANCODE_LANG2 },
    { "lang3", SDL_SCANCODE_LANG3 },
    { "lang4", SDL_SCANCODE_LANG4 },
    { "lang5", SDL_SCANCODE_LANG5 },
    { "lang6", SDL_SCANCODE_LANG6 },
    { "lang7", SDL_SCANCODE_LANG7 },
    { "lang8", SDL_SCANCODE_LANG8 },
    { "lang9", SDL_SCANCODE_LANG9 },

    { "alterase", SDL_SCANCODE_ALTERASE },
    { "sysreq", SDL_SCANCODE_SYSREQ },
    { "cancel", SDL_SCANCODE_CANCEL },
    { "clear", SDL_SCANCODE_CLEAR },
    { "prior", SDL_SCANCODE_PRIOR },
    { "return2", SDL_SCANCODE_RETURN2 },
    { "separator", SDL_SCANCODE_SEPARATOR },
    { "out", SDL_SCANCODE_OUT },
    { "oper", SDL_SCANCODE_OPER },
    { "clearagain", SDL_SCANCODE_CLEARAGAIN },
    { "crsel", SDL_SCANCODE_CRSEL },
    { "exsel", SDL_SCANCODE_EXSEL },

    { "kp_00", SDL_SCANCODE_KP_00 },
    { "kp_000", SDL_SCANCODE_KP_000 },
    { "thousandsseparator", SDL_SCANCODE_THOUSANDSSEPARATOR },
    { "decimalseparator", SDL_SCANCODE_DECIMALSEPARATOR },
    { "currencyunit", SDL_SCANCODE_CURRENCYUNIT },
    { "currencysubunit", SDL_SCANCODE_CURRENCYSUBUNIT },
    { "kp_leftparen", SDL_SCANCODE_KP_LEFTPAREN },
    { "kp_rightparen", SDL_SCANCODE_KP_RIGHTPAREN },
    { "kp_leftbrace", SDL_SCANCODE_KP_LEFTBRACE },
    { "kp_rightbrace", SDL_SCANCODE_KP_RIGHTBRACE },
    { "kp_tab", SDL_SCANCODE_KP_TAB },
    { "kp_backspace", SDL_SCANCODE_KP_BACKSPACE },
    { "kp_a", SDL_SCANCODE_KP_A },
    { "kp_b", SDL_SCANCODE_KP_B },
    { "kp_c", SDL_SCANCODE_KP_C },
    { "kp_d", SDL_SCANCODE_KP_D },
    { "kp_e", SDL_SCANCODE_KP_E },
    { "kp_f", SDL_SCANCODE_KP_F },
    { "kp_xor", SDL_SCANCODE_KP_XOR },
    { "kp_power", SDL_SCANCODE_KP_POWER },
    { "kp_percent", SDL_SCANCODE_KP_PERCENT },
    { "kp_less", SDL_SCANCODE_KP_LESS },
    { "kp_greater", SDL_SCANCODE_KP_GREATER },
    { "kp_ampersand", SDL_SCANCODE_KP_AMPERSAND },
    { "kp_dblampersand", SDL_SCANCODE_KP_DBLAMPERSAND },
    { "kp_verticalbar", SDL_SCANCODE_KP_VERTICALBAR },
    { "kp_dblverticalbar", SDL_SCANCODE_KP_DBLVERTICALBAR },
    { "kp_colon", SDL_SCANCODE_KP_COLON },
    { "kp_hash", SDL_SCANCODE_KP_HASH },
    { "kp_space", SDL_SCANCODE_KP_SPACE },
    { "kp_at", SDL_SCANCODE_KP_AT },
    { "kp_exclam", SDL_SCANCODE_KP_EXCLAM },
    { "kp_memstore", SDL_SCANCODE_KP_MEMSTORE },
    { "kp_memrecall", SDL_SCANCODE_KP_MEMRECALL },
    { "kp_memclear", SDL_SCANCODE_KP_MEMCLEAR },
    { "kp_memadd", SDL_SCANCODE_KP_MEMADD },
    { "kp_memsubtract", SDL_SCANCODE_KP_MEMSUBTRACT },
    { "kp_memmultiply", SDL_SCANCODE_KP_MEMMULTIPLY },
    { "kp_memdivide", SDL_SCANCODE_KP_MEMDIVIDE },
    { "kp_plusminus", SDL_SCANCODE_KP_PLUSMINUS },
    { "kp_clear", SDL_SCANCODE_KP_CLEAR },
    { "kp_clearentry", SDL_SCANCODE_KP_CLEARENTRY },
    { "kp_binary", SDL_SCANCODE_KP_BINARY },
    { "kp_octal", SDL_SCANCODE_KP_OCTAL },
    { "kp_decimal", SDL_SCANCODE_KP_DECIMAL },
    { "kp_hexadecimal", SDL_SCANCODE_KP_HEXADECIMAL },

    { "lctrl", SDL_SCANCODE_LCTRL },
    { "lshift", SDL_SCANCODE_LSHIFT },
    { "lalt", SDL_SCANCODE_LALT },
    { "lgui", SDL_SCANCODE_LGUI },
    { "rctrl", SDL_SCANCODE_RCTRL },
    { "rshift", SDL_SCANCODE_RSHIFT },
    { "ralt", SDL_SCANCODE_RALT },
    { "rgui", SDL_SCANCODE_RGUI },

    { "mode", SDL_SCANCODE_MODE },
    { "audionext", SDL_SCANCODE_AUDIONEXT },
    { "audioprev", SDL_SCANCODE_AUDIOPREV },
    { "audiostop", SDL_SCANCODE_AUDIOSTOP },
    { "audioplay", SDL_SCANCODE_AUDIOPLAY },
    { "audiomute", SDL_SCANCODE_AUDIOMUTE },
    { "mediaselect", SDL_SCANCODE_MEDIASELECT },
    { "www", SDL_SCANCODE_WWW },
    { "mail", SDL_SCANCODE_MAIL },
    { "calculator", SDL_SCANCODE_CALCULATOR },
    { "computer", SDL_SCANCODE_COMPUTER },
    { "ac_search", SDL_SCANCODE_AC_SEARCH },
    { "ac_home", SDL_SCANCODE_AC_HOME },
    { "ac_back", SDL_SCANCODE_AC_BACK },
    { "ac_forward", SDL_SCANCODE_AC_FORWARD },
    { "ac_stop", SDL_SCANCODE_AC_STOP },
    { "ac_refresh", SDL_SCANCODE_AC_REFRESH },
    { "ac_bookmarks", SDL_SCANCODE_AC_BOOKMARKS },

    { "brightnessdown", SDL_SCANCODE_BRIGHTNESSDOWN },
    { "brightnessup", SDL_SCANCODE_BRIGHTNESSUP },
    { "displayswitch", SDL_SCANCODE_DISPLAYSWITCH },

    { "kbdillumtoggle", SDL_SCANCODE_KBDILLUMTOGGLE },
    { "kbdillumdown", SDL_SCANCODE_KBDILLUMDOWN },
    { "kbdillumup", SDL_SCANCODE_KBDILLUMUP },
    { "eject", SDL_SCANCODE_EJECT },
    { "sleep", SDL_SCANCODE_SLEEP },

    { "app1", SDL_SCANCODE_APP1 },
    { "app2", SDL_SCANCODE_APP2 },
    // { "audiorewind", SDL_SCANCODE_AUDIOREWIND },
    // { "audiofastforward", SDL_SCANCODE_AUDIOFASTFORWARD },
};

}
