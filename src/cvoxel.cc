// Hacked-together C interface for writing to voxel worlds, for Python
// Ctypes. Not even going to think of Window's existence for now.

#include "voxels.hh"

#include <atomic>
#include <stdexcept>
#include <stdio.h>
#include <unordered_map>

static std::unique_ptr<myricube::VoxelWorld> p_world;

// Set the current world to modify by passing in the path to its
// world.myricube file.
extern "C" int myricube_select_world(const char* world_filename)
{
    try {
        p_world.reset(new myricube::VoxelWorld(world_filename));
        return 0;
    }
    catch (const std::runtime_error& e) {
        fprintf(stderr, "Opening '%s': %s\n", world_filename, e.what());
        return 1;
    }
    catch (...) {
        fprintf(stderr, "Unknown error opening '%s'\n", world_filename);
        return 1;
    }
}

// Set the voxel at the given coordinate in the current world.
extern "C" int myricube_set(int32_t x, int32_t y, int32_t z, uint32_t hex_color)
{
    if (p_world == nullptr) {
        fprintf(stderr, "No current world\n");
        return 1;
    }
    try {
        p_world->set(glm::ivec3(x,y,z), hex_color);
        return 0;
    }
    catch (...) {
        fprintf(stderr, "Unknown error in myricube_set\n");
        return 1;
    }
}

// Get the current voxel at the given coordinate in the current world.
// Warning: At the moment, this could create new chunk groups on disk
// (when we query a coordinate that's not in any current chunk group)
// so it's not a read-only operation.
extern "C" uint32_t myricube_get(int32_t x, int32_t y, int32_t z)
{
    if (p_world == nullptr) {
        fprintf(stderr, "No current world\n");
        return 0;
    }
    try {
        return (*p_world)(glm::ivec3(x,y,z));
    }
    catch (...) {
        fprintf(stderr, "Unknown error in myricube_get\n");
        return 0;
    }
}

