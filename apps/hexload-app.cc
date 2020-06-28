// Load the VoxelWorld (in myricube:hex format) named by the MYRICUBE_WORLD
// environment variable.

#include "app.hh"

#include <stdlib.h>

#include "hexload.hh"

namespace myricube {

void app_init(VoxelWorld& world, Window&)
{
    const char* filename = getenv("MYRICUBE_WORLD");
    if (filename == nullptr) {
        fprintf(stderr, "Missing MYRICUBE_WORLD environment variable.\n");
    }
    bool okay = read_hex(world, filename);
    if (!okay) {
        fprintf(stderr, "Load from %s failed.\n", filename);
    }
}

void app_update(VoxelWorld&)
{

}

} // end namespace
