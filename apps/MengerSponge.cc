#include "app.hh"

#include <thread>

namespace myricube {

template <int Levels>
class BaseMengerSponge : public App
{
    VoxelWorld world_;
    std::thread thread;

  public:
    BaseMengerSponge() : world_{}, thread(thread_loop, world_.get_handle())
    {
        thread.detach();
    }

    VoxelWorld& update(float) override
    {
        return world_;
    }

  private:
    static void thread_loop(WorldHandle world_handle)
    {
        VoxelWorld world(world_handle);
        int pow_3_level = 1;
        for (int i = 0; i < Levels; ++i) pow_3_level *= 3;
        f(world, Levels, pow_3_level, glm::ivec3(0, 0, 0), 0);
    }
    
    static void f(
        VoxelWorld& world,
        int level, int pow_3_level,
        glm::ivec3 offset,
        uint32_t color)
    {
        if (level == 0) {
            world.set(offset, color | 0x80000000);
            return;
        }
        int n = pow_3_level / 3;
        f(world, level-1, n, offset, 0x303030);
        f(world, level-1, n, offset + glm::ivec3(n,   0,   0), 0x803030);
        f(world, level-1, n, offset + glm::ivec3(2*n, 0,   0), 0xC03030);
        f(world, level-1, n, offset + glm::ivec3(2*n, n,   0), 0xC08030);
        f(world, level-1, n, offset + glm::ivec3(2*n, 2*n, 0), 0xC0C000);
        f(world, level-1, n, offset + glm::ivec3(n,   2*n, 0), 0x80C000);
        f(world, level-1, n, offset + glm::ivec3(0,   2*n, 0), 0x30C030);
        f(world, level-1, n, offset + glm::ivec3(0,   n,   0), 0x308030);

        f(world, level-1, n, offset + glm::ivec3(0,   0,   n), 0x303080);
        f(world, level-1, n, offset + glm::ivec3(2*n, 0,   n), 0xC03080);
        f(world, level-1, n, offset + glm::ivec3(2*n, 2*n, n), 0xC0C080);
        f(world, level-1, n, offset + glm::ivec3(0,   2*n, n), 0x30C080);

        f(world, level-1, n, offset + glm::ivec3(0,   0,   2*n), 0x3030C0);
        f(world, level-1, n, offset + glm::ivec3(n,   0,   2*n), 0x8030C0);
        f(world, level-1, n, offset + glm::ivec3(2*n, 0,   2*n), 0xC030C0);
        f(world, level-1, n, offset + glm::ivec3(2*n, n,   2*n), 0xC080C0);
        f(world, level-1, n, offset + glm::ivec3(2*n, 2*n, 2*n), 0xC0C0C0);
        f(world, level-1, n, offset + glm::ivec3(n,   2*n, 2*n), 0x80C0C0);
        f(world, level-1, n, offset + glm::ivec3(0,   2*n, 2*n), 0x30C0C0);
        f(world, level-1, n, offset + glm::ivec3(0,   n,   2*n), 0x3080C0);
    };
};

using MengerSponge = BaseMengerSponge<6>;

MYRICUBE_ADD_APP(MengerSponge)

} // end namespace
