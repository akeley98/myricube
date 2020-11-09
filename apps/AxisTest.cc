// Test axis code
#include "app.hh"

namespace myricube {

class AxisTest : public App
{
    VoxelWorld world;
  public:
    AxisTest()
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
            
            world.set(glm::ivec3(i, 1, 0), red);
            world.set(glm::ivec3(-i, 1, 0), red);
            world.set(glm::ivec3(0, i+1, 0), green);
            world.set(glm::ivec3(0, 1-i, 0), green);
            world.set(glm::ivec3(0, 1, i), blue);
            world.set(glm::ivec3(0, 1, -i), blue);
        }
    }
    
    VoxelWorld& update(float) override
    {
        return world;
    }
};

MYRICUBE_ADD_APP(AxisTest)

}
