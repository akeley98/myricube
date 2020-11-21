// 3D random walk generator.
#include "app.hh"

#include <random>

namespace myricube {

template <bool InMemory>
class RandomWalkBase : public App
{
    VoxelWorld world {
        InMemory
        ? add_in_memory_prefix("RandomWalk/world.myricube")
        : expand_filename("RandomWalk/world.myricube") };

    std::mt19937 rng{19980724};

    void add_random_walk()
    {
        uint8_t blue = rng() >> 24;
        uint8_t red = rng() >> 24;
        uint8_t green_base = 40 + (rng() >> 26);
        int x = 0, y = 0, z = 0;
        for (int i = 0; i < 888000; ++i) {
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
            world.set(glm::ivec3(x, y, z), voxel);
        }
    }

  public:
    RandomWalkBase()
    {
        add_random_walk();
    }

    VoxelWorld& update(float) override
    {
        return world;
    }

    void add_key_targets(Window& window) override
    {
        KeyTarget target;
        target.down = [this] (KeyArg arg) -> bool
        {
            if (arg.repeat) return false;
            this->add_random_walk();
            return true;
        };
        window.add_key_target("add_random_walk", target);

        KeyTarget target100;
        target100.down = [this] (KeyArg arg) -> bool
        {
            if (arg.repeat) return false;
            for (int i = 0; i < 100; ++i) this->add_random_walk();
            return true;
        };
        window.add_key_target("add_random_walk_100", target100);
    }
};

using RandomWalk = RandomWalkBase<false>;
using RandomWalkMem = RandomWalkBase<true>;

MYRICUBE_ADD_APP(RandomWalk)
MYRICUBE_ADD_APP(RandomWalkMem)

} // end namespace
