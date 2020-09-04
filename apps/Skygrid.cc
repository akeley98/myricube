// Fill in voxels where x/y/z coords are all multiples of 4.  This is
// probably a worst-case for raycasting -- lots of empty space to cast
// through.

#include "app.hh"

#include <random>

namespace myricube {

template <int Radius, int Spacing>
class BaseSkygrid : public App
{
    VoxelWorld world;

  public:
    BaseSkygrid()
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

        for (int z = -Radius; z <= Radius; z += Spacing) {
            for (int y = -Radius; y <= Radius; y += Spacing) {
                for (int x = -Radius; x <= Radius; x += Spacing) {
                    world.set(glm::ivec3(x,y,z), random_voxel());
                }
            }
        }
    }

    VoxelWorld& update(float) override
    {
        return world;
    }
};

class Skygrid : public BaseSkygrid<360, 4> {};
MYRICUBE_ADD_APP(Skygrid)

class SparseGrid : public BaseSkygrid<600, 32> {};
MYRICUBE_ADD_APP(SparseGrid)

}
