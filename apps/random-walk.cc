// 3D random walk generator.
#include "app.hh"

#include <random>

namespace myricube {

VoxelWorld* p_world;
std::mt19937 rng{19980724};

void add_random_walk()
{
    uint8_t blue = rng() >> 24;
    uint8_t red = rng() >> 24;
    uint8_t green_base = 40 + (rng() >> 26);
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
}

void app_init(VoxelWorld& world, Window& window)
{
    p_world = &world;
    KeyTarget target;
    target.down = [] (KeyArg arg) -> bool
    {
        if (arg.repeat) return false;
        add_random_walk();
        return true;
    };
    window.add_key_target("add_random_walk", target);
    add_random_walk();
}

void app_update(VoxelWorld& world)
{
    // I don't anticipate the world being relocated but to prevent trouble...
    p_world = &world;
}

} // end namespace
