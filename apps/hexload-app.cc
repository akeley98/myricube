// Load the VoxelWorld (in myricube:hex format) named by the myricube_world
// environment variable.

#include "app.hh"

#include <stdlib.h>

#include "hexload.hh"

namespace myricube {

class Hexload : public App
{
    VoxelWorld world;

  public:
    Hexload()
    {
        const char* filename = getenv("myricube_world");
        if (filename == nullptr) {
            panic("Missing myricube_world environment variable.");
        }
        bool okay = read_hex(world, filename);
        if (!okay) {
            fprintf(stderr, "Load from %s failed.\n", filename);
        }
    }

    VoxelWorld& update(float) override
    {
        return world;
    }
};

MYRICUBE_ADD_APP(Hexload)

} // end namespace
