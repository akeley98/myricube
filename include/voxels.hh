// CPU-side implementations of voxel grid data structures: chunks,
// chunk groups, and entire "infinite" voxel worlds stored on disk or
// in-memory and accessed via mmap. All files are accessed by a
// filename_string (8-bit on Unix, 16-bit on Windows); file names that
// start with the in_memory_prefix are (temporary) in-memory files.
//
// My rough plan is to have voxel worlds be stored in a single
// directory.  Each chunk group will be have its exact binary
// representation stored on disk in a file of the directory labeled
// with the group coordinate.
//
// Additionally, each chunk group in the "infinite" world (whether it
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
#include <vector>

namespace myricube {

// This is probably the data type you care about.
class VoxelWorld;



// File name of the BinWorld file, and filename format for chunk group files.
const char world_filename[] = "world.myricube";

constexpr uint64_t chunk_group_base_magic_number = 587569177;

constexpr uint64_t renderer_mesh_dirty_flag = uint64_t(1) << 63;
constexpr uint64_t renderer_raycast_dirty_flag = uint64_t(1) << 62;

// I changed the voxel endianness (red_shift, green_shift, blue_shift)
// partway through development, so I need this bit in the magic number
// to detect files in the old endianness.
constexpr uint64_t new_endian_magic = uint64_t(1) << 63;



// /MEM:/ is the in_memory_prefix. Colons are not allowed in file
// names on Windows and almost certainly not in the name of a
// top-level directory on Unix.
static const filename_string in_memory_prefix = "/MEM:/";

inline bool starts_with_in_memory_prefix(const filename_string& arg)
{
    if (arg.size() < in_memory_prefix.size()) return false;
    size_t index = 0;
    for (auto c : in_memory_prefix) {
        if (arg[index++] != c) return false;
    }
    return true;
}

inline filename_string add_in_memory_prefix(const filename_string& arg)
{
    if (starts_with_in_memory_prefix(arg)) return arg;

    filename_string result = in_memory_prefix;
    result.append(arg);
    return result;
}



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
        uint64_t(ChunkSize) << 40 |
        new_endian_magic;

    const uint64_t magic_number = expected_magic;
    uint64_t reserved[510] = { 0 };

    // Set to ~uint64_t(0) _after_ making modifications, or call
    // mark_dirty(). You can also set it periodically while doing
    // modifications to speed-up the renderer's view of it.
    //
    // If this chunk group is newly created on disk, remember to set
    // the corresponding BinGroupBitfield bit as well.
    mutable std::atomic<uint64_t> dirty_flags = { 0 };

    void mark_dirty() const
    {
        dirty_flags.store(~uint64_t(0));
    }

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

    // Given the GROUP COORDINATE of a group within the section of
    // chunk groups (i.e. this function masks the upper bits away for
    // you), return whether this bitfield records said chunk group as
    // existing on disk.
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
// Note that there's no actual link to chunk groups here; they are
// associated with this world file by being stored in the same
// directory.
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
        arg->mark_dirty();
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
// file system for (somewhat) handling synchronization problems.  The
// main purpose of this class is to open chunk group files on disk (or
// in the in-memory file system) and return memory mappings of them as
// smart pointers.
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
    // within the world's directory, OR the name of the world
    // directory with a trailing slash.
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

    // Overloaded version of mut/view_chunk_group.
    void get_chunk_group_uptr(
        glm::ivec3 group_coord, UPtrMutChunkGroup& ptr)
    {
        ptr = mut_chunk_group(group_coord);
    }

    void get_chunk_group_uptr(
        glm::ivec3 group_coord, UPtrChunkGroup& ptr) const
    {
        ptr = view_chunk_group(group_coord);
    }

    const filename_string& get_directory_name() const
    {
        return directory_trailing_slash;
    }
};



// Object for efficiently viewing and modifying chunks in a voxel
// world. Holds a direct-mapped cache of pointers to memory-mapped
// chunk groups, minimizing the need to interact with the file system.
//
// The object itself is not synchronized, but multiple caches for the
// same world may be used concurrently (assuming no OS bugs).  This is
// still not exactly super efficient -- the cache itself is a huge 1
// megabyte-ish structure.
template <typename UPtr, bool CheckBitfield>
class WorldCache
{
    WorldHandle world;

    static constexpr int32_t modulus = 32;

    // First, we need to store the bitfields corresponding to the
    // chunk groups we tried to read or modify. Every time we
    // encounter a new (or evicted) chunk group, we look into the
    // bitfield_vector for a pointer to its corresponding mmap'd
    // bitfield. If not found, we have to add it.
    //
    // Note that nothing is ever removed from bitfield_vector; in
    // practice this doesn't matter as bitfields cover a huge amount
    // of space, so this simplification is worth it to prevent
    // worrying about dangling pointers.
    struct CanonicalGroupBitfield
    {
        glm::ivec3 canonical_group_coord;
        UPtrGroupBitfield bitfield_ptr;
    };
    std::vector<CanonicalGroupBitfield> bitfield_vector;

  public:
    struct Entry
    {
        // Index into bitfield_vector that gives us the bitfield for
        // this entry. Negative value indicates this Entry has not
        // yet been populated (i.e. is invalid).
        int32_t bitfield_index = -1;

        // Group coordinate of the chunk group this Entry is for.
        glm::ivec3 group_coord;

        // Pointer to chunk group on disk. Null if it doesn't exist.
        UPtr chunk_group_ptr;

        // Helper for determining if this entry is for the group with
        // the given group coordinate.
        bool is_for_group(glm::ivec3 group_coord_arg) const
        {
            return group_coord == group_coord_arg and bitfield_index >= 0;
        }
    };
    Entry cache[modulus][modulus][modulus];

    explicit WorldCache(WorldHandle world_) : world(std::move(world_)) { }
    ~WorldCache() = default;
    WorldCache(const WorldCache&) = delete;

  private:
    // Return the location to store the Entry for the named chunk group.
    Entry& get_entry_location(glm::ivec3 group_coord)
    {
        return cache[group_coord.z & (modulus-1)]
                    [group_coord.y & (modulus-1)]
                    [group_coord.x & (modulus-1)];
    }

  public:
    // Get the Entry for the named chunk group, reading the named
    // chunk group and store it in the cache if needed.  The returned
    // reference is valid as long as this->get_entry is not called
    // again with a different group coordinate. (This is the reason
    // this is not a const member).
    const Entry& get_entry(glm::ivec3 group_coord)
    {
        Entry& entry = get_entry_location(group_coord);

        // Case 1: Already found in cache. In this case we just need to
        // check if there's a null pointer, and if so, check the bitfield
        // (if required) to see if the chunk group now exists on disk.
        if (entry.is_for_group(group_coord)) {
            if (entry.chunk_group_ptr == nullptr) {
                bool lookup = true;
                if (CheckBitfield) {
                    lookup = bitfield_vector.at(entry.bitfield_index)
                            .bitfield_ptr->chunk_group_on_disk(group_coord);
                }
                if (lookup) {
                    world.get_chunk_group_uptr(group_coord,
                                               entry.chunk_group_ptr);
                }
            }
            return entry;
        }

        // Case 2: Not in cache; we have to load it in.
        // Construct the new entry separately for exception safety.
        Entry new_entry;
        new_entry.group_coord = group_coord;

        // Search for the bitfield for this chunk group.
        new_entry.bitfield_index = 0;
        glm::ivec3 canonical =
            BinGroupBitfield::canonical_group_coord(group_coord);
        for (const auto& n : bitfield_vector) {
            if (n.canonical_group_coord == canonical) goto found;
            new_entry.bitfield_index++;
        }
        // If we reach (no goto) this point, have to load the
        // bitfield and APPEND to the bitfield vector (so
        // bitfield_index is correct).
        bitfield_vector.push_back(
            { canonical, world.bitfield_for_chunk_group(canonical) } );

      found:
        // Now check the bitfield (if required) to see if the chunk
        // group is on disk and if so load it in.
        bool lookup = true;
        if (CheckBitfield) {
            lookup = bitfield_vector.at(new_entry.bitfield_index)
                    .bitfield_ptr->chunk_group_on_disk(group_coord);
        }
        if (lookup) {
            world.get_chunk_group_uptr(group_coord, new_entry.chunk_group_ptr);
        }

        // No exceptions thrown: can safely overwrite evicted entry.
        entry = std::move(new_entry);
        return entry;
    }

    WorldHandle get_handle() const
    {
        return world;
    }
};

// Cache specialized for modifying voxel worlds: cache pointers to mutable,
// and skip bitfield checks (so nonexistent chunks will be created).
using MutWorldCache = WorldCache<UPtrMutChunkGroup, false>;

// Cache specialized for reading voxel worlds: cache pointers to const,
// and check bitfields to lower overhead of accessing nonexistent groups.
using ViewWorldCache = WorldCache<UPtrChunkGroup, true>;



// Friendly interface to a voxel world. Feed in the file name of the
// world file you want to open. By default, the "file name" is that of
// the temporary shared in-memory voxel world (not saved after program
// shutdown).
class VoxelWorld
{
    mutable MutWorldCache world_cache;
    static const inline filename_string in_memory_world_filename =
        filename_concat_c_str(in_memory_prefix, world_filename);

  public:
    explicit VoxelWorld(WorldHandle handle) :
        world_cache(handle) { }

    explicit VoxelWorld(filename_string filename=in_memory_world_filename) :
        VoxelWorld(WorldHandle(std::move(filename))) { }

    // Return voxel at the given coordinate.
    uint32_t operator() (glm::ivec3 c) const
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
        ptr->mark_dirty();
        assert(ptr != nullptr);
        ptr->set(c, voxel);
    }

    // For every voxel in the axis-aligned box with corners corner0
    // and corner1 (inclusive bounds), run the given callback on it.
    // Callback takes 2 args: a pointer to uint32_t (the voxel value)
    // and the coordinate of the voxel. The dirty flags for every
    // chunk group thus visited is set automatically.
    template <typename F>
    void map(glm::ivec3 corner0, glm::ivec3 corner1, F& f)
    {
        glm::ivec3 global_low = glm::min(corner0, corner1);
        glm::ivec3 global_high = glm::max(corner0, corner1);

        // Run f for every voxel in the given chunk that overlaps with
        // the caller-defined box. chunk_lower_left is the world
        // coordinate of the lower-left-most voxel in the chunk.
        auto map_chunk = [&] (BinChunk& chunk, glm::ivec3 chunk_lower_left)
        {
            glm::ivec3 local_low = glm::max(
                glm::ivec3(0), global_low - chunk_lower_left);
            glm::ivec3 local_high = glm::min(
                glm::ivec3(chunk_size-1), global_high - chunk_lower_left);

            for (int32_t z = local_low.z; z <= local_high.z; ++z) {
                for (int32_t y = local_low.y; y <= local_high.y; ++y) {
                    for (int32_t x = local_low.x; x <= local_high.x; ++x) {
                        glm::ivec3 coord = chunk_lower_left + glm::ivec3(x,y,z);
                        f(&chunk.voxel_array[z][y][x], coord);
                    }
                }
            }
        };

        // Given the group coordinate of a chunk group, map that chunk
        // and run f on the voxels within that chunk that overlap with
        // the caller-specified box. Also mark as dirty.
        auto map_group = [&] (glm::ivec3 group_coord)
        {
            BinChunkGroup* group =
                world_cache.get_entry(group_coord).chunk_group_ptr.get();

            for (int z = 0; z < edge_chunks; ++z) {
                for (int y = 0; y < edge_chunks; ++y) {
                    for (int x = 0; x < edge_chunks; ++x) {
                        BinChunk& chunk = group->chunk_array[z][y][x];
                        glm::ivec3 chunk_lower_left = group_coord * group_size
                                              + glm::ivec3(x,y,z) * chunk_size;
                        map_chunk(chunk, chunk_lower_left);
                    }
                }
            }

            group->mark_dirty();
        };

        // Map over the chunk groups that overlap with the caller's box.
        // Floor divide to get lower bound.
        int32_t xL = global_low.x >> group_shift;
        int32_t yL = global_low.y >> group_shift;
        int32_t zL = global_low.z >> group_shift;
        // Ceil divide to get upper bound.
        int32_t xH = (global_high.x + group_size - 1) >> group_shift;
        int32_t yH = (global_high.y + group_size - 1) >> group_shift;
        int32_t zH = (global_high.z + group_size - 1) >> group_shift;

        for (int32_t x = xL; x <= xH; ++x) {
            for (int32_t y = yL; y <= yH; ++y) {
                for (int32_t z = zL; z <= zH; ++z) {
                    map_group(glm::ivec3(x, y, z));
                }
            }
        }
    }

    // These next two functions are potentially dangerous, but the
    // costs of the friendly functions above are high enough that I
    // don't want to force usage of them.
    MutWorldCache& borrow_cache()
    {
        return world_cache;
    }

    WorldHandle get_handle()
    {
        return world_cache.get_handle();
    }
};

} // end namespace myricube

#endif /* !MYRICUBE_VOXELS_HH_ */
