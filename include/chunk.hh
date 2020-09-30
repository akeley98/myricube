// Semi-temporary compatibility class. Interface to a mmap'd
// scratch/temporary voxel world that emulates my older in-memory
// voxel world interface. chunk.hh name is historical.

#ifndef MYRICUBE_CHUNK_HH_
#define MYRICUBE_CHUNK_HH_

#include "myricube.hh"

#include <assert.h>
#include <atomic>

#include "voxels.hh"

namespace myricube {

class Renderer;

// Compatibility handle to the global on-disk "scratch" voxel world.
class VoxelWorld
{
    friend class Renderer;

    WorldHandle handle = { expand_filename("scratch/world.myricube") };
    MutWorldCache world_cache = { handle };

  public:
    // Return voxel at the given coordinate.
    uint32_t operator() (glm::ivec3 c)
    {
        glm::ivec3 group_coord = to_group_coord(c);
        BinChunkGroup* ptr =
            world_cache.get_entry(group_coord).chunk_group_ptr.get();
        return ptr ? (*ptr)(c) : 0;
    }

    // Set the voxel at the given coordinate.
    void set(glm::ivec3 c, uint32_t voxel)
    {
        glm::ivec3 group_coord = to_group_coord(c);
        BinChunkGroup* ptr =
            world_cache.get_entry(group_coord).chunk_group_ptr.get();
        ptr->dirty_flags.store(~uint64_t(0));
        assert(ptr != nullptr);
        ptr->set(c, voxel);
    }
};

} // end namespace
#endif /* !MYRICUBE_CHUNK_HH_ */
