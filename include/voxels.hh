// CPU-side implementations of voxel grid data structures: chunks,
// chunk groups, and entire "infinite" voxel worlds stored on disk and
// accessed via mmap.
//
// My rough plan is to have voxel worlds be stored in a single
// directory.  Each chunk group will be have its exact binary
// representation stored on disk in a file of the directory labeled
// with the group coordinate.
//
// Additionally, each chunk group is the "infinite" world (whether it
// actually exists or not) is assigned a bit in a group bitfield that
// tells whether the chunk group actually exists on disk or not --
// this is to save time on repeatedly looking up files for chunk
// groups that do not actually exist yet.
//
// NOTE: Some structs are templatized in case I want to make changes to the
// chunk / group sizes later and need to convert world files as a result.

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

constexpr uint64_t chunk_group_base_magic_number = 587569177;

constexpr uint64_t renderer_mesh_dirty_flag = uint64_t(1) << 63;
constexpr uint64_t renderer_raycast_dirty_flag = uint64_t(1) << 62;



// Voxels are really represented as 32-bit packed ints now. This is
// needed for compatibility.
struct Voxel
{
    bool visible = false;
    uint8_t red = 0, green = 0, blue = 0;
    operator bool() const { return visible; }
    operator uint32_t() const
    {
        return uint32_t(blue) << blue_shift
             | uint32_t(green) << green_shift
             | uint32_t(red) << red_shift
             | (visible ? visible_bit : 0);
    }

    Voxel() = default;
    Voxel(const Voxel&) = default;
    Voxel& operator= (const Voxel&) = default;

    Voxel(uint8_t red_, uint8_t green_, uint8_t blue_)
    {
        visible = true;
        red = red_;
        green = green_;
        blue = blue_;
    }

    Voxel(uint32_t packed)
    {
        visible = (packed & visible_bit) != 0;
        red   = (packed >> red_shift)   & 255u;
        green = (packed >> green_shift) & 255u;
        blue  = (packed >> blue_shift)  & 255u;
    }
};



// Chunk of voxels as it appears in binary on-disk.
template <size_t ChunkSize>
struct BinChunkT
{
    // Voxels (packed 32-bit color) in this chunk, in [z][y][x] order.
    uint32_t voxel_array[ChunkSize][ChunkSize][ChunkSize] = { };

    // Given the world or residue or in-chunk coordinate of a voxel
    // within this chunk (masking makes all those coordinates the
    // same), return said voxel.
    uint32_t operator() (glm::ivec3 coord) const
    {
        return voxel_array[coord.z & (ChunkSize-1)]
                          [coord.y & (ChunkSize-1)]
                          [coord.x & (ChunkSize-1)];
    }

    // Set the voxel with the given world/residue/chunk coordinate.
    void set(glm::ivec3 coord, uint32_t voxel)
    {
        voxel_array[coord.z & (ChunkSize-1)]
                   [coord.y & (ChunkSize-1)]
                   [coord.x & (ChunkSize-1)] = voxel;
    }
};

using BinChunk = BinChunkT<chunk_size>;



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

    // Set to ~uint64_t(0) _after_ making modifications.  You can also
    // set it periodically while doing modifications to speed-up the
    // renderer's view of it.
    //
    // If this chunk group is newly created on disk, remember to set
    // the corresponding BinGroupBitfield bit as well.
    mutable std::atomic<uint64_t> dirty_flags = { 0 };

    // Chunks within this chunk group, in [z][y][x] order.
    BinChunkT<ChunkSize> chunk_array[EdgeChunks][EdgeChunks][EdgeChunks];

    // Given the world or residue coordinate of a voxel in this chunk
    // group (again masking equalizes all those systems), return said
    // voxels' value.
    uint32_t operator() (glm::ivec3 coord) const
    {
        auto group_mask = EdgeChunks * ChunkSize - 1;
        auto residue_x = coord.x & group_mask;
        auto residue_y = coord.y & group_mask;
        auto residue_z = coord.z & group_mask;
        return chunk_array[residue_z / ChunkSize]
                          [residue_y / ChunkSize]
                          [residue_x / ChunkSize] (coord);
    }

    // Like above but set the voxel's value.
    void set(glm::ivec3 coord, uint32_t voxel)
    {
        auto group_mask = EdgeChunks * ChunkSize - 1;
        auto residue_x = coord.x & group_mask;
        auto residue_y = coord.y & group_mask;
        auto residue_z = coord.z & group_mask;
        return chunk_array[residue_z / ChunkSize]
                          [residue_y / ChunkSize]
                          [residue_x / ChunkSize].set(coord, voxel);
    }
};

using BinChunkGroup = BinChunkGroupT<edge_chunks, chunk_size>;

static_assert(sizeof(BinChunkGroup) ==
    4096 + sizeof(uint32_t) * (group_size * group_size * group_size));



// File on-disk for recording whether chunk groups exist on disk or not.
// Each 64 x 64 x 64 section of chunk groups share one BinGroupBitfield.
struct BinGroupBitfield
{
    static constexpr int32_t modulus = 64;

    // Read next function for format.
    std::atomic<uint64_t> x_bitfield_zy[modulus][modulus] = {};

    // Given the GROUP COORDINATE of a group within the section of chunk
    // groups (i.e. upper bits masked away), return whether this bitfield
    // records said chunk group as existing on disk.
    bool chunk_group_on_disk(glm::ivec3 group_coord)
    {
        auto x = group_coord.x & (modulus - 1);
        auto y = group_coord.y & (modulus - 1);
        auto z = group_coord.z & (modulus - 1);
        return (x_bitfield_zy[z][y].load() >> x) & 1u;
    }

    // Set the bit corresponding to the above.
    void set_chunk_group_on_disk(glm::ivec3 group_coord)
    {
        auto x = group_coord.x & (modulus - 1);
        auto y = group_coord.y & (modulus - 1);
        auto z = group_coord.z & (modulus - 1);
        x_bitfield_zy[z][y] |= uint64_t(1) << x;
    }

    // Return the group coordinate of the lower-left chunk group
    // covered by this group bitfield. Gives this bitfield a unique
    // identifier.
    static glm::ivec3 canonical_group_coord(glm::ivec3 group_coord)
    {
        group_coord.x &= ~(modulus - 1);
        group_coord.y &= ~(modulus - 1);
        group_coord.z &= ~(modulus - 1);
        return group_coord;
    }
};



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



// Similar deal for the chunk group bitfield.
void unmap_bin_group_bitfield(BinGroupBitfield*);
struct GroupBitfieldDeleter {
    void operator () (BinGroupBitfield* arg)
    {
        unmap_bin_group_bitfield(arg);
    }
};
using UPtrGroupBitfield =
    std::unique_ptr<BinGroupBitfield, GroupBitfieldDeleter>;



// Handle for a voxel world (directory storing chunk groups).
// Multiple handles to the same world may exist; we rely on the OS
// file system for (somewhat) handling synchronization problems.
class WorldHandle
{
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
    uint64_t get_last_nz_atomic() const
    {
        assert(bin_world != nullptr);
        return bin_world->atomic_counter.load();
    }

    // Return the result of the "next" call to inc_nz_atomic.
    uint64_t get_next_nz_atomic() const
    {
        assert(bin_world != nullptr);
        return 1u + bin_world->atomic_counter.load();
    }

    // These next functions are rather expensive as they involve file
    // system access. Typically, use the WorldMutator class instead,
    // which caches results.

    // Return a pointer suitable for modifying the memory-mapped chunk
    // group with the given group coordinate. The chunk group is
    // created on disk if it does not exist.
    UPtrMutChunkGroup mut_chunk_group(glm::ivec3);

    // Return a pointer for viewing the chunk group with the given
    // group coordinate. Returns nullptr if there is no such chunk
    // group.
    UPtrChunkGroup view_chunk_group(glm::ivec3) const;

    // Return the entire group bitfield (not an individual bit!) that
    // holds the bit for the given chunk group.
    UPtrGroupBitfield bitfield_for_chunk_group(glm::ivec3) const;
};

} // end namespace myricube

#endif /* !MYRICUBE_VOXELS_HH_ */
