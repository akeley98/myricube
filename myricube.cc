// See header file for overview.

#include "myricube.hh"

#include <assert.h>
#include <stdio.h>

#include <random>
#include <stdexcept>
#include <string>
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
    return data_directory + in;
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
}

void bind_keys(Window& window)
{
    window.bind_keycode(SDL_SCANCODE_LEFT, "pop_old_camera");
    window.bind_keycode(-8, "pop_old_camera");
    window.bind_keycode(SDL_SCANCODE_RIGHT, "pop_future_camera");
    window.bind_keycode(-9, "pop_future_camera");

    window.bind_keycode(SDL_SCANCODE_U, "forward");
    window.bind_keycode(SDL_SCANCODE_SPACE, "backward");
    window.bind_keycode(SDL_SCANCODE_P, "leftward");
    window.bind_keycode(SDL_SCANCODE_A, "rightward");
    window.bind_keycode(SDL_SCANCODE_LCTRL, "upward");
    window.bind_keycode(SDL_SCANCODE_LALT, "downward");

    window.bind_keycode(SDL_SCANCODE_O, "sprint");
    window.bind_keycode(SDL_SCANCODE_I, "speed_up");
    window.bind_keycode(SDL_SCANCODE_COMMA, "slow_down");

    window.bind_keycode(-3, "look_around");
    window.bind_keycode(-4, "vertical_scroll");
    window.bind_keycode(-5, "vertical_scroll");
    window.bind_keycode(-6, "horizontal_scroll");
    window.bind_keycode(-7, "horizontal_scroll");

    window.bind_keycode(SDL_SCANCODE_K, "do_it");
    window.bind_keycode(SDL_SCANCODE_K, "add_random_walk");
    window.bind_keycode(SDL_SCANCODE_Z, "pause");
    window.bind_keycode(SDL_SCANCODE_B, "toggle_chunk_debug");
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
        void draw_skybox(glm::mat4, glm::mat4);
        draw_skybox(camera.get_residue_view(), camera.get_projection());
        render_world_mesh_step(world, camera);
        render_world_raycast_step(world, camera);
        window.set_title("Myricube "
                         + std::to_string(window.get_fps()) + " FPS");
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
