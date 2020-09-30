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

    // Last-used pointer for reading from a chunk group, and its group coord.
    // Used to minimize looking up files.
    struct PositionedReadChunkGroup
    {
        glm::ivec3 group_coord;
        UPtrChunkGroup ptr;
    };
    PositionedReadChunkGroup cached_read;

    // Similar to above, but for writing.
    struct PositionedWriteChunkGroup
    {
        glm::ivec3 group_coord;
        UPtrMutChunkGroup ptr;
    };
    PositionedWriteChunkGroup cached_write;

    WorldHandle handle = { expand_filename("scratch/world.myricube") };

  public:
    // Meaningless placeholder for now; all VoxelWorld refer to the
    // same world.
    uint64_t id() const
    {
        return 1;
    }

    // Return voxel at the given coordinate.
    uint32_t operator() (glm::ivec3 c)
    {
        glm::ivec3 group_coord = to_group_coord(c);
        if (group_coord != cached_read.group_coord
        or cached_read.ptr == nullptr)
        {
            cached_read.ptr = handle.view_chunk_group(group_coord);
            cached_read.group_coord = group_coord;
        }
        return cached_read.ptr ? (*cached_read.ptr)(group_coord) : 0;
    }

    // Set the voxel at the given coordinate.
    void set(glm::ivec3 c, uint32_t voxel)
    {
        glm::ivec3 group_coord = to_group_coord(c);
        if (group_coord != cached_write.group_coord
        or cached_write.ptr == nullptr)
        {
            cached_write.ptr = handle.mut_chunk_group(group_coord);
            cached_write.group_coord = group_coord;
        }
        return cached_write.ptr->set(c, voxel);
    }
};

} // end namespace
#endif /* !MYRICUBE_CHUNK_HH_ */
