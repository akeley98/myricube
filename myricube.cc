// See header file for overview.

#include "myricube.hh"

#include <assert.h>
#include <stdio.h>

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "camera.hh"
#include "chunk.hh"
#include "renderer.hh"
#include "window.hh"

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

    Window window(on_window_resize);
    window.set_title("Myricube");

    VoxelWorld world;
    Camera camera;
    Voxel red(255, 0, 0);

    for (int i = 2; i < 20; i += 3) {
        world.set(glm::ivec3(i,0,0), red);
        world.set(glm::ivec3(-i,0,0), red);
        world.set(glm::ivec3(0,i,0), red);
        world.set(glm::ivec3(0,-i,0), red);
        world.set(glm::ivec3(0,0,i), red);
        world.set(glm::ivec3(0,0,-i), red);
    }

    gl_first_time_setup();
    while (window.update_swap_buffers(5)) {
        gl_clear();
        render_world_mesh_step(world, camera);
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
