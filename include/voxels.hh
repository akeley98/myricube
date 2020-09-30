// CPU-side implementations of voxel grid data structures: chunks,
// chunk groups, and entire "infinite" voxel worlds stored on disk and
// accessed via mmap.
//
// My rough plan is to have voxel worlds be stored in a single directory.
// Each chunk group will be have its exact binary representation stored
// on disk in a file of the directory labeled with the group coordinate.

#ifndef MYRICUBE_VOXELS_HH_
#define MYRICUBE_VOXELS_HH_

#include "myricube.hh"

#include <atomic>
#include <memory>
#include <stdint.h>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <string.h>

namespace myricube {

// File name of the BinWorld file, and filename format for chunk group files.
const char world_filename[] = "world.myricube";

// Windows uses UTF-16, sigh. (wstring is 16-bit on windows). TODO test.
#if defined(__WIN32__) || defined(__WIN32) || defined(WIN32)
#define MYRICUBE_WINDOWS 1
using filename_string = std::wstring;
#endif

// Meanwhile the free world uses UTF-8 like civilized individuals.
#if defined(__linux__) || defined(__unix__)
using filename_string = std::string;
#endif

using filename_char = filename_string::value_type;

// Dumb function for appending an ascii C string to a filename.
// Needed to deal with platform differences.
inline filename_string filename_concat_c_str(
    filename_string filename, const char* c_str)
{
    auto len = strlen(c_str);
    filename.reserve(filename.size() + len + 1);

    for (size_t i = 0; i < len; ++i) {
        char c = c_str[i];
        assert(0 < c and c <= 127); // Must be ASCII.
        filename.push_back(static_cast<filename_char>(c));
    }
    return filename;
}

constexpr uint64_t chunk_group_base_magic_number = 587569177;

// Chunk of voxels as it appears in binary on-disk.
template <size_t ChunkSize>
struct BinChunkT
{
    // Voxels (packed 32-bit color) in this chunk, in [z][y][x] order.
    uint32_t voxel_array[ChunkSize][ChunkSize][ChunkSize] = { };
};

// Chunk group as it appears in binary on-disk.
template <size_t EdgeChunks, size_t ChunkSize>
struct BinChunkGroupT
{
    static constexpr uint64_t expected_magic =
        chunk_group_base_magic_number |
        uint64_t(EdgeChunks) << 32 |
        uint64_t(ChunkSize) << 40;

    const uint64_t magic_number = expected_magic;
    uint64_t reserved[510] = { 0 };

    // Set to UINT64_MAX _after_ making modifications.
    // You can also set it periodically while doing modifications
    // to speed-up the renderer's view of it.
    std::atomic<uint64_t> dirty_flags = { 0 };

    // Chunks within this chunk group, in [z][y][x] order.
    BinChunkT<ChunkSize> chunk_array[EdgeChunks][EdgeChunks][EdgeChunks];
};

using BinChunkGroup = BinChunkGroupT<edge_chunks, chunk_size>;

static_assert(sizeof(BinChunkGroup) ==
    4096 + sizeof(uint32_t) * (group_size * group_size * group_size));

// File on-disk for representing the whole voxel world.
//
// Note that there's no actual link to chunk groups here; they
template <size_t EdgeChunks, size_t ChunkSize>
struct BinWorldT
{
    const char shebang[80] =
        "#!/usr/bin/env myricube\n";

    static constexpr uint64_t expected_magic =
        BinChunkGroupT<EdgeChunks, ChunkSize>::expected_magic;

    const uint64_t magic_number = expected_magic;

    // Should only ever increase. This is forever, so don't increase
    // it by more than needed (although with 64 bits this should not
    // overflow in practice). WorldHandle methods recommended in
    // preference to manual modification.
    std::atomic<uint64_t> atomic_counter = { 0x1'0000'0000 };

    uint64_t reserved[500] = { 0 };
};

using BinWorld = BinWorldT<edge_chunks, chunk_size>;

static_assert(sizeof(BinWorld) == 4096);

// Deleter for BinChunkGroup returned by WorldHandle (to be declared).
// You can rely on mutable chunk groups having their dirty flag set
// upon deletion.
void unmap_mut_bin_chunk_group(BinChunkGroup*);
void unmap_bin_chunk_group(const BinChunkGroup*);
struct MutChunkGroupDeleter {
    void operator () (BinChunkGroup* arg)
    {
        arg->dirty_flags.store(~uint64_t(0));
        unmap_mut_bin_chunk_group(arg);
    }
};

struct ChunkGroupDeleter {
    void operator () (const BinChunkGroup* arg)
    {
        unmap_bin_chunk_group(arg);
    }
};

// Non-const and const smart pointers for BinChunkGroup. They are
// intentionally not interchangable (different deleters).
using UPtrMutChunkGroup = std::unique_ptr<BinChunkGroup, MutChunkGroupDeleter>;
using UPtrChunkGroup = std::unique_ptr<const BinChunkGroup, ChunkGroupDeleter>;

// Handle for a voxel world (directory storing chunk groups).
// Multiple handles to the same world may exist; we rely on the OS
// file system for (somewhat) handling synchronization problems.
class WorldHandle
{
    friend void unmap_mut_bin_chunk_group(BinChunkGroup*);
    friend void unmap_bin_chunk_group(const BinChunkGroup*);

    // Directory that the world is stored in, with trailing slash
    // so filenames can be appended directly.
    filename_string directory_trailing_slash;

    // Pointer to memory-mapped BinWorld file. shared_ptr custom
    // deleter (passed at runtime elsewhere unlike unique_ptr for some
    // reason) will handle closing the memory mapping.
    std::shared_ptr<BinWorld> bin_world;

  public:
    // Construct world handle for world stored in the given directory.
    // The actual filename passed is the path of the BinWorld file
    // within the world's directory.
    WorldHandle(const filename_string& world_filename);

    ~WorldHandle() = default;
    WorldHandle(WorldHandle&&) = default;
    WorldHandle(const WorldHandle&) = default;
    WorldHandle& operator=(WorldHandle&&) = default;
    WorldHandle& operator=(const WorldHandle&) = default;

    // Return a pointer suitable for modifying the memory-mapped chunk
    // group with the given group coordinate. The chunk group is
    // created on disk if it does not exist.
    UPtrMutChunkGroup mut_chunk_group(glm::ivec3);

    // Return a pointer for viewing the chunk group with the given
    // group coordinate. Returns nullptr if there is no such chunk
    // group.
    UPtrChunkGroup view_chunk_group(glm::ivec3);

    // Each world has its own 64-bit monotonic-increasing atomic
    // counter.  (Assuming no-one goes behind our back and modifies it
    // on disk inappropriately).

    // Increment and return the new atomic counter value. It's
    // going to be nonzero (for a very long time), so 0 can be
    // used as a sentinel value.
    uint64_t inc_nz_atomic()
    {
        assert(bin_world != nullptr);
        return ++(bin_world->atomic_counter);
    }

    // Return the result of the "previous" call to inc_nz_atomic.
    uint64_t get_last_nz_atomic()
    {
        assert(bin_world != nullptr);
        return bin_world->atomic_counter.load();
    }

    // Return the result of the "next" call to inc_nz_atomic.
    uint64_t get_next_nz_atomic()
    {
        assert(bin_world != nullptr);
        return 1u + bin_world->atomic_counter.load();
    }
};

} // end namespace myricube

#endif /* !MYRICUBE_VOXELS_HH_ */
