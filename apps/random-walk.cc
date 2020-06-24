// 3D random walk generator.
#include "app.hh"

#include <random>

namespace myricube {

VoxelWorld* p_world;
std::mt19937 rng;

void app_init(VoxelWorld& world, Window& window)
{
    p_world = &world;
    KeyTarget target;
    target.down = [] (KeyArg arg) -> bool
    {
        if (arg.repeat) return false;
        uint8_t blue = 127 + (rng() >> 25);
        uint8_t red = 127 + (rng() >> 25);
        uint8_t green_base = rng() >> 25;
        int x = 0, y = 0, z = 0;
        for (int i = 0; i < 200000; ++i) {
            switch (rng() % 6) {
                case 0: x++; break;
                case 1: y++; break;
                case 2: z++; break;
                case 3: x--; break;
                case 4: y--; break;
                case 5: z--; break;
            }
            uint8_t green = uint8_t((rng() >> 25) + green_base);
            Voxel voxel(red, green, blue);
            p_world->set(glm::ivec3(x, y, z), voxel);
        }
        return true;
    };
    window.add_key_target("add_random_walk", target);
}

void app_update(VoxelWorld& world)
{
    // I don't anticipate the world being relocated but to prevent trouble...
    p_world = &world;
}

} // end namespace
