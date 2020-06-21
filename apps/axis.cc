// Test axis code
#include "myricube.hh"
#include "app.hh"

namespace myricube {

void app_init(VoxelWorld& world, Window&)
{
    for (int i = 0; i <= 128; i += 4) {
        uint8_t n = 0;
        if (i % 16 == 0) n = 128;
        if (i % 64 == 0) n = 255;
        Voxel red(255, n, n);
        Voxel green(n, 255, n);
        Voxel blue(n, n, 255);
        
        world.set(glm::ivec3(i, 0, 0), red);
        world.set(glm::ivec3(-i, 0, 0), red);
        world.set(glm::ivec3(0, i, 0), green);
        world.set(glm::ivec3(0, -i, 0), green);
        world.set(glm::ivec3(0, 0, i), blue);
        world.set(glm::ivec3(0, 0, -i), blue);
    }
}

void app_update(VoxelWorld&)
{

}

}
