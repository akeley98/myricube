// See header file for overview.

#include "myricube.hh"

#include <assert.h>
#include <stdio.h>

#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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

std::mt19937 rng;

void add_random_column(VoxelWorld& world)
{
    auto x = uint8_t(rng() % 128);
    auto y = uint8_t(rng() % 128);
    auto z = uint8_t(rng() % 128);
    auto r = uint8_t(rng() % 256);
    auto b = uint8_t(rng() % 256);
    for (int i = 0; i < 80; ++i) {
        world.set(glm::ivec3(x, y+i, z), Voxel(r,i*3,b));
    }
}

void add_key_targets(Window& window, Camera& camera, VoxelWorld& world)
{
    static float speed = 8.0f;
    static float sprint_mod = 1.0f;

    KeyTarget forward, backward, leftward, rightward, upward, downward;
    forward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera.frenet_move(0, 0, +arg.dt * speed * sprint_mod);
        return true;
    };
    backward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera.frenet_move(0, 0, -arg.dt * speed * sprint_mod);
        return true;
    };
    leftward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera.frenet_move(-arg.dt * speed * sprint_mod, 0, 0);
        return true;
    };
    rightward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera.frenet_move(+arg.dt * speed * sprint_mod, 0, 0);
        return true;
    };
    upward.per_frame = [&] (KeyArg arg) -> bool
    {
        camera.frenet_move(0, +arg.dt * speed * sprint_mod, 0);
        return true;
    };
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

    KeyTarget add_random_block;
    add_random_block.down = [&] (KeyArg arg) -> bool
    {
        add_random_column(world);
        return true;
    };
    window.add_key_target("add_random_block", add_random_block);
}

void bind_keys(Window& window)
{
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

    window.bind_keycode(SDL_SCANCODE_K, "add_random_block");
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
    add_key_targets(window, camera, world);
    bind_keys(window);

    Voxel red(180, 0, 100);
    Voxel green(0, 130, 0);
    Voxel blue(140, 200, 255);
    Voxel yellow(255, 255, 0);

    for (int i = 0; i < 100; i += 3) {
        auto red_or_yellow = i % 30 == 0 ? yellow : red;
        auto green_or_yellow = i % 30 == 0 ? yellow : green;
        auto blue_or_yellow = i % 30 == 0 ? yellow : blue;
        world.set(glm::ivec3(i,0,0), red_or_yellow);
        world.set(glm::ivec3(-i,0,0), yellow);
        world.set(glm::ivec3(0,i,0), green_or_yellow);
        world.set(glm::ivec3(0,-i,0), yellow);
        world.set(glm::ivec3(0,0,i), blue_or_yellow);
        world.set(glm::ivec3(0,0,-i), yellow);
    }

    gl_first_time_setup();
    while (window.update_swap_buffers(5)) {
        gl_clear();
        void draw_skybox(glm::mat4, glm::mat4);
        draw_skybox(camera.get_residue_view(), camera.get_projection());
        render_world_mesh_step(world, camera);
        add_random_column(world);
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
