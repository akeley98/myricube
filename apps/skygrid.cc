// Fill in voxels where x/y/z coords are all multiples of 4.  This is
// probably a worst-case for raycasting -- lots of empty space to cast
// through.

#include "app.hh"

#include <random>

namespace myricube {

// Called once at the start; also an opportunity to add key targets.
void app_init(VoxelWorld& world, Window&)
{
    std::mt19937 rng;
    
    auto random_voxel = [&rng] () -> Voxel
    {
        auto n = rng();
        uint8_t red = (n >> 8) & 255;
        uint8_t green = (n >> 16) & 255;
        uint8_t blue = (n >> 24) & 255;
        return Voxel(red, green, blue);
    };
    
    constexpr int radius = 512;
    for (int z = -radius; z <= radius; z += 4) {
        for (int y = -radius; y <= radius; y += 4) {
            for (int x = -radius; x <= radius; x += 4) {
                world.set(glm::ivec3(x,y,z), random_voxel());
            }
        }
    }
}

// Called per frame.
void app_update(VoxelWorld&)
{

}

}
