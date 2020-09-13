// Implementation of the hybrid voxel renderer. I don't see how I can
// make much of this exception safe, so I wrap everything in a
// noexcept at the end. As usual with OpenGL, none of this is
// thread-safe either.

#include "myricube.hh"

#include <algorithm>
#include <memory>
#include <stdio.h>
#include <typeinfo>
#include <utility>
#include <vector>

#include "camera.hh"
#include "chunk.hh"
#include "glad/glad.h"
#include "renderer.hh"
#include "shaders.hh"

namespace myricube {

// Hacky debug variables.
bool chunk_debug = false;
bool evict_stats_debug = false;
bool disable_zcull_sort = false;

// Extract color from voxel and pack as 32-bit integer.
inline uint32_t to_packed_color(Voxel v)
{
    return uint32_t(v.blue) << blue_shift
         | uint32_t(v.green) << green_shift
         | uint32_t(v.red) << red_shift
         | (v.visible ? visible_bit : 0);
}

// Write out the voxel as a texel in the given GL format and type (if
// supported).
template <GLenum Format, GLenum Type>
inline void write_voxel_texel(Voxel, void*)
{
    static_assert(Format != Format, "Add support for format.");
}

template<> inline void write_voxel_texel<GL_RGBA, GL_UNSIGNED_INT_8_8_8_8>
    (Voxel v, void* texel)
{
    *static_cast<uint32_t*>(texel) =
        uint32_t(v.red) << 24 |
        uint32_t(v.green) << 16 |
        uint32_t(v.blue) << 8 |
        uint32_t(v.visible ? 255 : 0);
}

// My new plan for the GPU voxel mesh: each visible voxel will be
// represented by ONE vertex in a VBO; instanced rendering will
// transform this vertex into an actual cube. This vertex encodes, in
// a packed bitfield format,
//
// The residue coordinates of the voxel
//
// The color of the voxel (8-bits each for R, G, B)
//
// Which faces of the voxel are visible (+/- x, y, z faces).
// (This depends on whether the neighboring voxels are filled or not)
struct MeshVoxelVertex
{
    uint32_t packed_residue_face_bits = 0xFFFFFFFF;
    uint32_t packed_color = 0xFFFFFFFF;

    MeshVoxelVertex() = default;

    // Given a _visible_ voxel and its residue coordinate, return the
    // VBO vertex carrying this information.
    //
    // This constructor defines the packed layout.
    MeshVoxelVertex(Voxel v, uint8_t x, uint8_t y, uint8_t z)
    {
        static_assert(group_size <= 255,
                      "group too big for 8-bit unsigned coordinates.");
        packed_residue_face_bits = uint32_t(x) << x_shift
                                 | uint32_t(y) << y_shift
                                 | uint32_t(z) << z_shift;
        packed_color = to_packed_color(v);
        assert(packed_color & visible_bit);
    }
};

// The maximum number of MeshVoxelVerts needed for one chunk.  For now
// it's just the total number of voxels per chunk. I'm sure there's a
// lower possible bound but for now I'll be conservative (even though
// I'm desperate for GPU memory).
constexpr size_t chunk_max_verts = chunk_size * chunk_size * chunk_size;

// CPU-side copy of the "mesh" (i.e. list of visible voxels) for a
// single chunk.
struct ChunkMesh
{
    MeshVoxelVertex verts[chunk_max_verts];
    size_t vert_count = 0;
};

// Handle for GPU resources for the mesh of one chunk group.  The
// CPU-side state of this MeshEntry should always accurately describe
// the state of the copy of this data on the GPU and GL. (i.e. if you
// modify the data within, it is your responsibility to update the GPU
// state at the same time).
struct MeshEntry
{
    // World ID and group coordinate of the chunk group this mesh is
    // for. This can be used to tell if this entry is correct or needs
    // to be replaced in MeshStore.
    uint64_t world_id = 0;
    glm::ivec3 group_coord;

    // mesh_array[z][y][x] is the mesh for ChunkGroup::chunk_array[z][y][x].
    // This is copied into the VBO.
    ChunkMesh mesh_array[edge_chunks][edge_chunks][edge_chunks];

    // OpenGL name for the VBO used to store the meshes of this chunk group.
    // 0 when not yet allocated.
    GLuint vbo_name = 0;

    // Size in bytes of the vbo's data store on the GPU.
    static constexpr GLsizeiptr vbo_bytes = sizeof mesh_array;

    // Return the offset (in number of MeshVoxelVertexes, not bytes)
    // into the VBO where the data from mesh_array[z][y][x] is copied
    // into.
    static unsigned vert_offset(unsigned x, unsigned y, unsigned z)
    {
        assert(x < edge_chunks and y < edge_chunks and z < edge_chunks);
        unsigned chunk_idx = x + y*edge_chunks + z*edge_chunks*edge_chunks;
        return chunk_idx * chunk_max_verts;
    }

    // Same as above, but return offset as count of bytes.
    static GLsizeiptr byte_offset(unsigned x, unsigned y, unsigned z)
    {
        auto vert_sz = GLsizeiptr(sizeof(ChunkMesh::verts[0]));
        GLsizeiptr off = vert_offset(x, y, z) * vert_sz;
        assert(size_t(off) < vbo_bytes);
        return off;
    }

    MeshEntry() = default;

    // Get some RAII going.
    ~MeshEntry()
    {
        glDeleteBuffers(1, &vbo_name);
        vbo_name = 0;
        world_id = 0;
    }

    MeshEntry(MeshEntry&& other)
    {
        swap(*this, other);
    }

    friend void swap(MeshEntry& left, MeshEntry& right)
    {
        using std::swap;
        swap(left.world_id, right.world_id);
        swap(left.group_coord, right.group_coord);
        swap(left.mesh_array, right.mesh_array);
        swap(left.vbo_name, right.vbo_name);
    }
};

// For memory efficiency, the AABB of a chunk is stored in packed
// format on the GPU.
struct PackedAABB
{
    uint32_t packed_low;
    uint32_t packed_high;

    static_assert(group_size <= 255,
                  "group too big for 8-bit unsigned coordinates.");

    PackedAABB() = default;

    PackedAABB(glm::ivec3 aabb_low, glm::ivec3 aabb_high)
    {
        auto x = aabb_low.x;
        auto y = aabb_low.y;
        auto z = aabb_low.z;
        packed_low = uint32_t(x) << x_shift
                   | uint32_t(y) << y_shift
                   | uint32_t(z) << z_shift;
        x = aabb_high.x;
        y = aabb_high.y;
        z = aabb_high.z;
        packed_high = uint32_t(x) << x_shift
                    | uint32_t(y) << y_shift
                    | uint32_t(z) << z_shift;
    }

    uint8_t low_x() const
    {
        return uint8_t((packed_low >> x_shift) & 255);
    }
    uint8_t low_y() const
    {
        return uint8_t((packed_low >> y_shift) & 255);
    }
    uint8_t low_z() const
    {
        return uint8_t((packed_low >> z_shift) & 255);
    }
    uint8_t high_x() const
    {
        return uint8_t((packed_high >> x_shift) & 255);
    }
    uint8_t high_y() const
    {
        return uint8_t((packed_high >> y_shift) & 255);
    }
    uint8_t high_z() const
    {
        return uint8_t((packed_high >> z_shift) & 255);
    }
};

// All voxels within a chunk group share a single 3D voxel array on
// the GPU (stored as a 3D texture), as well as one VBO used to store
// the array of AABB for the chunks within the group. This is the
// handle for the GPU data needed to raycast one chunk group.
struct RaycastEntry
{
    // World ID and group coordinate of the chunk group this texture
    // and AABB array is for. This can be used to tell if this entry
    // is correct or needs to be replaced in RaycastStore.
    uint64_t world_id = 0;
    glm::ivec3 group_coord;

    // OpenGL name of the VBO storing the aabb_array. 0 when not yet
    // allocated.
    GLuint vbo_name = 0;

    // OpenGL name for the 3D voxels texture. 0 when not yet allocated.
    GLuint texture_name = 0;

    // Used for my staging buffer texture upload scheme later.
    bool needs_memory_barrier = false;
    bool should_draw = false;

    RaycastEntry() = default;

    ~RaycastEntry()
    {
        GLuint buffers[1] = { vbo_name };
        glDeleteBuffers(1, buffers);
        vbo_name = 0;
        glDeleteTextures(1, &texture_name);
        texture_name = 0;
    }

    RaycastEntry(RaycastEntry&& other)
    {
        swap(*this, other);
    }

    friend void swap(RaycastEntry& left, RaycastEntry& right)
    {
        using std::swap;
        swap(left.world_id, right.world_id);
        swap(left.group_coord, right.group_coord);
        swap(left.vbo_name, right.vbo_name);
        swap(left.texture_name, right.texture_name);
        swap(left.needs_memory_barrier, right.needs_memory_barrier);
        swap(left.should_draw, right.should_draw);
    }
};

// Struct for requesting a chunk group's GPU data from a MeshStore /
// RaycastStore entry. Some stats are also returned.
struct StoreRequest
{
    // Input: allow the store to return a nullptr if servicing this
    // request would be slow somehow.
    bool may_fail = false;

    // Input: allow the store to overwrite data in order to service
    // this request. If set, may_fail must also be set.
    bool read_only = false;

    // Output: whether the requested chunk group was already in the
    // store (not necessarily up-to-date).
    bool was_in_store;
};

// Cache of MeshEntry/RaycastEntry objects. The idea is to store entry
// structs for chunk groups near the camera. They are stored in a
// wraparound fashion (i.e.  using modular arithmetic on group
// coordinates) so that as the camera moves, the new groups being
// rendered overwrites the old groups no longer being rendered.
template <class Entry, uint32_t DefaultModulus, uint32_t Assoc>
class BaseStore
{
    struct CacheSet
    {
        uint64_t last_access[Assoc] = {0};
        Entry slots[Assoc];
    };
    // Monotonic increasing counter for least-recently accessed cache
    // eviction algorithm.
    uint64_t access_counter = 0;

    // Used for tracking how many cache entries were read for this frame.
    uint64_t frame_begin_access_counter = 0;

    // Cache is (modulus x modulus x modulus x Assoc) in size.
    uint32_t modulus = DefaultModulus;

    // Conceptually, a 3D (modulus x modulus x modulus) array of cache
    // sets (CacheSet itself comprises the Assoc dimension).
    //
    // Declared after modulus defensively.
    std::unique_ptr<CacheSet[]> cache_sets {
        new CacheSet[DefaultModulus * DefaultModulus * DefaultModulus] };

    // Fully-associative victim cache. Evicted Entry objects are
    // swapped into here.
    struct VictimCacheSlot
    {
        uint64_t last_access = 0;
        Entry entry;
    };
    std::vector<VictimCacheSlot> victim_cache;

    // Return a pointer to the location that the Entry for the chunk
    // group with the given group coordinate and world id should be.
    // The p_valid bool tells us whether the returned Entry matches
    // the one being searched for.
    Entry* cached_location(
        bool* p_valid,
        glm::ivec3 group_coord,
        uint64_t world_id,
        bool may_fail=false)
    {
        uint32_t x = uint32_t(group_coord.x ^ 0x8000'0000) % modulus;
        uint32_t y = uint32_t(group_coord.y ^ 0x8000'0000) % modulus;
        uint32_t z = uint32_t(group_coord.z ^ 0x8000'0000) % modulus;
        CacheSet& cache_set = cache_sets[
            z * modulus * modulus + y * modulus + x];

        // Try to search for a valid Entry for the requested world &
        // group coordinate. Mark it as accessed if so.
        for (unsigned i = 0; i < Assoc; ++i) {
            Entry* entry = &cache_set.slots[i];
            if (entry->world_id == world_id
            and entry->group_coord == group_coord) {
                cache_set.last_access[i] = ++access_counter;
                *p_valid = true;
                return entry;
            }
        }

        // Not found at this point; need to choose the least-recently
        // used entry to evict. This is much colder code than above.
        // However, skip all this if we don't actually need it
        // (may_fail is true).
        if (may_fail) {
            *p_valid = false;
            return nullptr;
        }
        uint64_t min_access = cache_set.last_access[0];
        unsigned evict_idx = 0;
        for (unsigned i = 1; i < Assoc; ++i) {
            auto this_last_access = cache_set.last_access[i];
            if (this_last_access < min_access) {
                evict_idx = i;
                min_access = this_last_access;
            }
            if (evict_stats_debug) {
                fprintf(stderr, "%u %lu\n", i, this_last_access);
            }
        }
        Entry* evict_entry = &cache_set.slots[evict_idx];
        auto old_access_counter = cache_set.last_access[evict_idx];
        cache_set.last_access[evict_idx] = ++access_counter;
        if (evict_stats_debug) {
            fprintf(stderr, "Evict [%u][%u][%u] %u\n\n", x, y, z, evict_idx);
            if (old_access_counter > frame_begin_access_counter) {
                fprintf(stderr, "\x1b[35m\x1b[1mCACHE THRASHING WARNING\n"
                    "ChunkGroup (%i, %i, %i) and (%i, %i, %i) tried to use\n"
                    "the same cache slot in one frame. If the victim cache\n"
                    "is not big enough, this may lead to flickering\n"
                    "(both will try to use the same memory in one frame, and\n"
                    "this bypasses the per-frame synchronization.)\x1b[0m\n",
                    evict_entry->group_coord.x,
                    evict_entry->group_coord.y,
                    evict_entry->group_coord.z,
                    group_coord.x, group_coord.y, group_coord.z);
            }
        }

        // Now we need to check the victim cache (if it exists),
        // either for the Entry being searched for (which will be
        // swapped back out of the victim cache), or for the oldest
        // victim cache entry to overwrite.
        bool found_in_victim_cache = false;
        if (!victim_cache.empty()) {
            VictimCacheSlot* victim_slot = &victim_cache[0];
            min_access = victim_slot->last_access;

            for (VictimCacheSlot& slot : victim_cache) {
                if (slot.entry.world_id == world_id
                and slot.entry.group_coord == group_coord) {
                    victim_slot = &slot;
                    found_in_victim_cache = true;
                    break;
                }
                else if (slot.last_access < min_access) {
                    victim_slot = &slot;
                    min_access = slot.last_access;
                }
            }
            using std::swap;
            swap(victim_slot->entry, *evict_entry);
            victim_slot->last_access = ++access_counter;
        }

        *p_valid = found_in_victim_cache;
        return evict_entry;
    }

  public:
    uint64_t eviction_count = 0;
    // If possible, return a valid, updated Entry for the given chunk
    // group extracted from the given world. See StoreRequest for
    // control parameters.
    //
    // This requires a template parameter providing static functions:
    //
    // replace(PositionedChunkGroup, VoxelWorld, Entry*)
    //
    //     Fill the Entry with data for this new chunk group,
    //     initializing OpenGL resources if needed.
    //
    // bool update(PositionedChunkGroup, VoxelWorld, Entry*, bool read_only)
    //
    //     Assuming replace was already called with identical
    //     arguments, above, just re-upload any changed data, unless
    //     read_only is true. Return true iff the Entry is
    //     up-to-date (i.e. return false only when read_only was
    //     true, but some changed data needed to be re-uploaded).
    template <typename EntryFiller>
    Entry* request(PositionedChunkGroup& pcg,
                   VoxelWorld& world,
                   StoreRequest* request)
    {
        assert(!request->read_only or request->may_fail);

        glm::ivec3 gc = group_coord(pcg);
        bool valid;
        Entry* entry =
            cached_location(&valid, gc, world.id(), request->may_fail);
        request->was_in_store = valid;

        // If the Entry in the array already corresponds to the given
        // chunk group; just update it if we're allowed to (deal with
        // dirty chunks) and return.
        if (valid) {
            bool success =
                EntryFiller::update(pcg, world, entry, request->read_only);
            assert(success or request->read_only);
            return success ? entry : nullptr;
        }

        // Otherwise, either fail if allowed, or evict and replace a
        // cache entry with data for the new chunk group.
        if (request->may_fail) return nullptr;

        EntryFiller::replace(pcg, world, entry);
        assert(entry->world_id == world.id() and entry->group_coord == gc);
        ++eviction_count;
        return entry;
    }

    // Thingie bolted-on after the fact to stretch my BaseStore class
    // far past its original design goals. This is just returns a
    // location for an entry with the given world/group_coord, but no
    // update/initialization functions are run (except filling in the
    // world_id and group_coord). StoreRequest handled as above.
    //
    // The raycast and mesh storage strategies are different enough
    // now that I probably should stop trying to share the common code
    // (or really factor out just the common LRU algorithm).
    Entry* request_no_update(glm::ivec3 gc,
                             uint64_t world_id,
                             StoreRequest* request)
    {
        assert(!request->read_only or request->may_fail);

        bool valid;
        Entry* entry =
            cached_location(&valid, gc, world_id, request->may_fail);
        request->was_in_store = valid;

        // If the Entry in the array already corresponds to the given
        // chunk group; just update it if we're allowed to (deal with
        // dirty chunks) and return.
        if (valid) {
            return entry;
        }

        // Otherwise, either fail if allowed, or evict and replace a
        // cache entry with data for the new chunk group.
        if (request->may_fail) return nullptr;

        entry->world_id = world_id;
        entry->group_coord = gc;
        ++eviction_count;
        return entry;
    }

    // Shrink/grow the cache, moving every cache entry (belonging to
    // the specified world) to its correct location in the new cache.
    // Entries for other worlds are skipped over.
    //
    // I put a noexcept because I don't really know how to properly
    // handle exceptions in the middle of OpenGL code, not because
    // this is guaranteed foolproof.
    void set_modulus(uint32_t new_modulus, uint64_t world_id) noexcept
    {
        if (new_modulus == modulus) return;
        assert(int(new_modulus) > 0);

        const uint32_t old_modulus = modulus;
        modulus = new_modulus;

        // Swap the new and old caches.
        std::unique_ptr<CacheSet[]> old_cache_sets(
            new CacheSet[new_modulus * new_modulus * new_modulus]);
        using std::swap;
        swap(old_cache_sets, cache_sets);

        fprintf(stderr, "Set modulus to %i\n", int(modulus));

        // Put the old entries in their correct place.
        if (world_id != 0) {
            const size_t end = old_modulus * old_modulus * old_modulus;
            for (size_t i = 0; i < end; ++i) {
                CacheSet& old_set = old_cache_sets[i];
                for (size_t a = 0; a < Assoc; ++a) {
                    Entry& old_entry = old_set.slots[a];
                    if (old_entry.world_id != world_id) continue;

                    bool valid;
                    Entry* new_location = cached_location(
                        &valid, old_entry.group_coord, world_id);
                    swap(*new_location, old_entry);
                    // Note: May want to preserve the old last_accessed time.
                    // For now, each cached_location updates the access time,
                    // so earlier repositioned cache entries are arbitrarily
                    // considered lower priority later.
                }
            }
        }
    }

    // Set the size of the victim cache. Same disclaimer about
    // noexcept applies.
    void resize_victim_cache(uint32_t cache_size) noexcept
    {
        // I've abandoned the victim cache for now (but for reasons
        // I'll keep to myself I'm going to keep it here instead of
        // relying on git history). Remove this assert and re-test if
        // I want it again.
        assert(cache_size == 0);
        victim_cache.resize(cache_size);
    }

    // Performance debug code. Used to measure how many cache entries
    // were read in this frame.
    void stats_begin_frame()
    {
        eviction_count = 0;
        frame_begin_access_counter = ++access_counter;
    }

    void print_stats_end_frame()
    {
        for (unsigned z = 0; z < modulus; ++z) {
            for (unsigned y = 0; y < modulus; ++y) {
                for (unsigned x = 0; x < modulus; ++x) {
                    auto& sets =
                        cache_sets[z * modulus * modulus + y * modulus + x];
                    unsigned count = 0;
                    for (size_t i = 0; i < Assoc; ++i) {
                        auto access = sets.last_access[i];
                        count += (access > frame_begin_access_counter);
                    }
                    if (count == Assoc) {
                        fprintf(stderr, "\x1b[1m\x1b[31m");
                    }
                    fprintf(stderr, "[%u][%u][%u] %u/%u  ",
                        z, y, x, count, Assoc);
                    for (unsigned i = 0; i < count; ++i) fprintf(stderr, "*");
                    fprintf(stderr, "\x1b[0m\n");
                }
            }
        }
        unsigned count = 0;
        for (const VictimCacheSlot& slot : victim_cache) {
            count += (slot.last_access > frame_begin_access_counter);
        }
        if (count == victim_cache.size()) {
            fprintf(stderr, "\x1b[1m\x1b[31m");
        }
        fprintf(stderr, "Victim Cache %u/%u\x1b[0m\n",
            count, unsigned(victim_cache.size()));
        fprintf(stderr, "%i evictions\n", int(eviction_count));
    }
};

// To be explained later (maybe).
struct StagingBuffer
{
    // World ID and group coordinate of the chunk group being staged
    // here. world_id = 0 indicates nothing is being staged here.
    uint64_t world_id = 0;
    glm::ivec3 group_coord;

    // aabb_array[z][y][x] is the minimal AABB containing the visible
    // voxels of ChunkGroup::chunk_array[z][y][x].
    //
    // If you change this, be aware that sizeof is used on this to
    // know the amount of bytes to allocate for the GPU's VBO.
    struct AABBs
    {
        PackedAABB aabb_array[edge_chunks][edge_chunks][edge_chunks];
    };

    // OpenGL name of the VBO storing the aabb_array. 0 when not yet
    // allocated.
    GLuint vbo_name = 0;

    // OpenGL name for the 3D voxels texture. 0 when not yet allocated.
    GLuint texture_name = 0;

    // OpenGL name of the staging SSBO.
    GLuint ssbo_name = 0;

    // Persistent mapped SSBO pointer.
    ChunkGroupVoxels* mapped_ssbo = nullptr;

    friend void swap(StagingBuffer& left, StagingBuffer& right)
    {
        using std::swap;
        swap(left.world_id, right.world_id);
        swap(left.group_coord, right.group_coord);
        swap(left.vbo_name, right.vbo_name);
        swap(left.texture_name, right.texture_name);
        swap(left.ssbo_name, right.ssbo_name);
        swap(left.mapped_ssbo, right.mapped_ssbo);
    }

    StagingBuffer() = default;

    StagingBuffer(StagingBuffer&& other)
    {
        swap(*this, other);
    }

    StagingBuffer& operator= (StagingBuffer&& other)
    {
        swap(*this, other);
        return *this;
    }

    ~StagingBuffer()
    {
        // Un-needed 0 check helps the C++ compiler.
        if (vbo_name != 0) {
            glDeleteBuffers(1, &vbo_name);
        }
        if (texture_name != 0) {
            glDeleteTextures(1, &texture_name);
        }
        if (ssbo_name != 0) {
            glDeleteBuffers(1, &ssbo_name);
            mapped_ssbo = nullptr;
        }
    }

    void make_ready()
    {
        if (vbo_name == 0) {
            glCreateBuffers(1, &vbo_name);
            auto sz = sizeof(AABBs);
            auto flags = GL_DYNAMIC_STORAGE_BIT;
            glNamedBufferStorage(vbo_name, sz, nullptr, flags);
            PANIC_IF_GL_ERROR;
        }

        if (texture_name == 0) {
            glCreateTextures(GL_TEXTURE_3D, 1, &texture_name);

            glBindTexture(GL_TEXTURE_3D, texture_name);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
            PANIC_IF_GL_ERROR;
            glTexStorage3D(
                GL_TEXTURE_3D, 1, GL_RGBA8, group_size, group_size, group_size);
            PANIC_IF_GL_ERROR;
        }

        if (ssbo_name == 0) {
            glCreateBuffers(1, &ssbo_name);
            auto sz = sizeof(uint32_t) * group_size * group_size * group_size;
            auto flags = GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT;
            glNamedBufferStorage(ssbo_name, sz, nullptr, flags);
            flags = GL_MAP_PERSISTENT_BIT | GL_MAP_INVALIDATE_BUFFER_BIT
                  | GL_MAP_FLUSH_EXPLICIT_BIT | GL_MAP_WRITE_BIT;
            mapped_ssbo = static_cast<ChunkGroupVoxels*>(glMapNamedBufferRange(
                ssbo_name, 0, sz, flags));
            PANIC_IF_GL_ERROR;
        }
    }
};

struct QueueEntry
{
    uint64_t world_id;
    glm::ivec3 group_coord;
};

class MeshStore : public BaseStore<MeshEntry, 2, 8> { };

class RaycastStore : public BaseStore<RaycastEntry, 4, 12>
{
  public:
    std::vector<StagingBuffer> write_staging_buffers, read_staging_buffers;
    std::vector<QueueEntry> queue;
};

// AABB is drawn as a unit cube from (0,0,0) to (1,1,1), which is
// stretched and positioned to the right shape and position in space.
// The normals are "flat shaded" and needed for dealing with
// floor/ceil rounding errors. (Basically, I add a bit of the normal
// vector to each face in order to ensure the voxels fit comfortably
// within the AABB).
//
// [positions] [normal]
//
// This really should be an array of structs, but I copied this from
// some older project and I'm too scared to screw it up.
static const float unit_box_vertices[48] =
{
    0, 1, 1,   -1, 0, 0,    // 0 -x face provoking vertex
    0, 0, 1,   0, 0, 1,     // 1 +z face provoking vertex
    1, 0, 1,   1, 0, 0,     // 2 +x face provoking vertex
    1, 1, 1,   0, 0, 0,     // 3 unused as provoking vertex
    0, 1, 0,   0, 1, 0,     // 4 +y face provoking vertex
    0, 0, 0,   0, -1, 0,    // 5 -y face provoking vertex
    1, 0, 0,   0, 0, 0,     // 6 unused as provoking vertex
    1, 1, 0,   0, 0, -1,    // 7 -z face provoking vertex
};

static const GLushort unit_box_elements[36] = {
    6, 7, 2, 7, 3, 2,   // +x face
    4, 5, 0, 5, 1, 0,   // -x face
    0, 3, 4, 3, 7, 4,   // +y face
    6, 2, 5, 2, 1, 5,   // -y face
    2, 3, 1, 3, 0, 1,   // +z face
    6, 5, 7, 5, 4, 7,   // -z face
};


// Voxels rendered using the mesh renderer also need a unit cube, but
// this time, instead of a normal I need a bit indicating which face
// each vertex is part of (so I can discard hidden faces).
struct VoxelUnitBoxVertex
{
    int32_t face_bit;
    float x, y, z;
    float u, v;
};

static const VoxelUnitBoxVertex voxel_unit_box_vertices[24] =
{
    { neg_x_face_bit, 0, 1, 1, 1, 1 },
    { neg_x_face_bit, 0, 1, 0, 1, 0 },
    { neg_x_face_bit, 0, 0, 1, 0, 1 },
    { neg_x_face_bit, 0, 0, 0, 0, 0 },

    { pos_x_face_bit, 1, 0, 0, 0, 0 },
    { pos_x_face_bit, 1, 1, 0, 1, 0 },
    { pos_x_face_bit, 1, 0, 1, 0, 1 },
    { pos_x_face_bit, 1, 1, 1, 1, 1 },

    { neg_y_face_bit, 0, 0, 0, 0, 0 },
    { neg_y_face_bit, 1, 0, 1, 1, 1 },
    { neg_y_face_bit, 0, 0, 1, 0, 1 },
    { neg_y_face_bit, 1, 0, 0, 1, 0 },

    { pos_y_face_bit, 1, 1, 1, 1, 1 },
    { pos_y_face_bit, 1, 1, 0, 1, 0 },
    { pos_y_face_bit, 0, 1, 0, 0, 0 },
    { pos_y_face_bit, 0, 1, 1, 0, 1 },

    { neg_z_face_bit, 0, 0, 0, 0, 0 },
    { neg_z_face_bit, 0, 1, 0, 0, 1 },
    { neg_z_face_bit, 1, 0, 0, 1, 0 },
    { neg_z_face_bit, 1, 1, 0, 1, 1 },

    { pos_z_face_bit, 0, 1, 1, 0, 1 },
    { pos_z_face_bit, 1, 0, 1, 1, 0 },
    { pos_z_face_bit, 1, 1, 1, 1, 1 },
    { pos_z_face_bit, 0, 0, 1, 0, 0 },
};

static const GLushort voxel_unit_box_elements[36] =
{
    0, 1, 2,    1, 3, 2,
    4, 5, 6,    5, 7, 6,
    8, 9, 10,   8, 11, 9,
    12, 13, 14, 15, 12, 14,
    16, 17, 18, 19, 18, 17,
    20, 21, 22, 21, 20, 23,
};

// Renderer class. Instantiate it with the camera and world to render
// and use it once.
class Renderer
{
    Camera& camera;
    VoxelWorld& world;
  public:
    Renderer(VoxelWorld& world_, Camera& camera_) :
        camera(camera_), world(world_) { }

    // Culling strategy: There are two phases to culling chunks: group
    // culling and distance culling. In group culling, entire groups of
    // chunks are culled if they are outside the view frustum, or if they
    // contain no visible voxels. Then of the chunks in non-culled groups,
    // their distance to the eye is calculated, and the chunk is culled,
    // drawn with raycasting, or drawn as a mesh, depending on whether
    // it is out-of-range of, far, or near the camera.

    // Should this group be culled? If not, return through
    // *squared_dist the squared distance between the eye and the
    // nearest point in the group.
    bool cull_group(PositionedChunkGroup& pcg, float* squared_dist=nullptr)
    {
        if (group(pcg).total_visible == 0) return true;

        glm::mat4 vp = camera.get_residue_vp();
        glm::ivec3 eye_group;
        glm::vec3 eye_residue;
        camera.get_eye(&eye_group, &eye_residue);

        // Never cull the group the eye is in (this avoids
        // pathological cases e.g. 7 corners behind eye and 1 in front
        // and out of the frustum).
        if (group_coord(pcg) == eye_group) {
            if (squared_dist) *squared_dist = 0.0f;
            return false;
        }

        eye_residue = glm::floor(eye_residue); // to match decide_chunk.

        // Position of this chunk group relative to the group that the
        // eye is in, in voxel units.
        glm::vec3 low_corner(group_size * (group_coord(pcg) - eye_group));

        const glm::vec3 x_edge = glm::vec3(group_size, 0, 0);
        const glm::vec3 y_edge = glm::vec3(0, group_size, 0);
        const glm::vec3 z_edge = glm::vec3(0, 0, group_size);

        // Compute the clip space (?) coordinates of the 8 corners of
        // this chunk group, and find the min/max of the xyz coordinates.
        // Take the abs of w before dividing x and y (but not z) so that
        // the x/y bounding planes don't swap sides.
        bool first_time = true;
        glm::vec3 low(-1), high(-1);
        auto minmax_corner = [vp, &low, &high, &first_time] (glm::vec3 v)
        {
            auto vp_v = vp * glm::vec4(v, 1);
            float abs_w = glm::abs(vp_v.w);
            glm::vec3 clip_coord(vp_v.x/abs_w, vp_v.y/abs_w, vp_v.z/vp_v.w);

            if (first_time) {
                low = clip_coord;
                high = clip_coord;
                first_time = false;
            }
            else {
                low = glm::min(low, clip_coord);
                high = glm::max(high, clip_coord);
            }
        };

        minmax_corner(low_corner);
        minmax_corner(low_corner + x_edge);
        minmax_corner(low_corner + y_edge);
        minmax_corner(low_corner + z_edge);
        minmax_corner(low_corner + x_edge + y_edge);
        minmax_corner(low_corner + x_edge + z_edge);
        minmax_corner(low_corner + y_edge + z_edge);
        glm::vec3 high_corner = low_corner + x_edge + y_edge + z_edge;
        minmax_corner(high_corner);

        // Cull the chunk if it is entirely out-of-bounds in one
        // direction on the x/y/z axis (one-direction == we won't clip
        // the group if the corners are all out-of-bounds but some
        // visible portion of the group "straddles" the frustum).
        if (low.x < -1 && high.x < -1) return true;
        if (low.y < -1 && high.y < -1) return true;
        if (low.z < 0 && high.z < 0) return true;
        if (low.x > 1 && high.x > 1) return true;
        if (low.y > 1 && high.y > 1) return true;
        if (low.z > 1 && high.z > 1) return true;

        if (squared_dist) {
            glm::vec3 nearest_point = glm::clamp(
                eye_residue, low_corner, high_corner);
            glm::vec3 disp = nearest_point - eye_residue;
            *squared_dist = glm::dot(disp, disp);
        }
        return false;
    }

    static constexpr int cull = 1, draw_mesh = 2, draw_raycast = 3;
    // Distance culling step. Given the group coordinate and AABB of a
    // chunk, decide which of the above things we should do: cull it,
    // draw it with the mesh algorithm (near the eye), or draw it with
    // the raycast algorithm (far away).
    //
    // This function may be duplicated on the GPU, hence the
    // floor(eye) -- this prevents subtle disagreements due to
    // different FP representations.
    int decide_chunk(glm::ivec3 group_coord,
                     Chunk& chunk)
    {
        if (chunk.total_visible == 0) return cull;
        glm::ivec3 aabb_low, aabb_high;
        chunk.get_aabb(&aabb_low, &aabb_high);
        glm::ivec3 eye_group;
        glm::vec3 eye_residue;
        camera.get_eye(&eye_group, &eye_residue);
        auto far_plane = camera.get_far_plane();
        auto raycast_thresh = camera.get_raycast_threshold();

        // Compute the displacement between the eye and the nearest point
        // of the chunk's AABB (to the eye).
        glm::ivec3 floor_eye = glm::ivec3(glm::floor(eye_residue))
                             + group_size * (eye_group - group_coord);
        glm::ivec3 aabb_nearest = glm::clamp(floor_eye, aabb_low, aabb_high);
        auto disp = glm::vec3(aabb_nearest - floor_eye);

        float squared_dist = glm::dot(disp, disp);
        if (squared_dist > far_plane * far_plane) return cull;
        if (squared_dist < raycast_thresh * raycast_thresh) return draw_mesh;
        return draw_raycast;
    }

    // Given a chunk, fill in the mesh needed to render it. I also
    // need the residue coordinate of the lower-left of the chunk
    // (i.e.  the position of the chunk relative to the chunk group it
    // is in) because the mesh is defined to be in residue
    // coordinates.
    static void fill_mesh_verts(ChunkMesh& mesh,
                                const Chunk& chunk,
                                glm::ivec3 chunk_residue)
    {
        // Look up whether the voxel at the given coordinate
        // (relative to the lower-left of this chunk) is visible.
        // Act as if voxels outside the chunk are always invisible.
        auto visible_block = [&chunk] (glm::ivec3 coord) -> bool
        {
            if (coord.x < 0 or coord.x >= chunk_size
             or coord.y < 0 or coord.y >= chunk_size
             or coord.z < 0 or coord.z >= chunk_size) return false;
            // Note: the masking in Chunk::operator () won't mess up
            // this coord. (I'm kind of violating my own comment in
            // Chunk::operator() because coord won't actually be in
            // the chunk unless that chunk is at (0,0,0)).
            return chunk(coord).visible;
        };

        auto visit_voxel = [&mesh, &chunk, visible_block, chunk_residue]
        (glm::ivec3 coord)
        {
            Voxel v = chunk(coord);
            if (!v.visible) return;

            uint8_t x = uint8_t(coord.x) + chunk_residue.x;
            uint8_t y = uint8_t(coord.y) + chunk_residue.y;
            uint8_t z = uint8_t(coord.z) + chunk_residue.z;

            MeshVoxelVertex vert(v, x, y, z);

            // Check which of the six faces are visible.
            if (!visible_block(coord + glm::ivec3(-1, 0, 0))) {
                vert.packed_residue_face_bits |= neg_x_face_bit;
            }
            if (!visible_block(coord + glm::ivec3(1, 0, 0))) {
                vert.packed_residue_face_bits |= pos_x_face_bit;
            }
            if (!visible_block(coord + glm::ivec3(0, -1, 0))) {
                vert.packed_residue_face_bits |= neg_y_face_bit;
            }
            if (!visible_block(coord + glm::ivec3(0, 1, 0))) {
                vert.packed_residue_face_bits |= pos_y_face_bit;
            }
            if (!visible_block(coord + glm::ivec3(0, 0, -1))) {
                vert.packed_residue_face_bits |= neg_z_face_bit;
            }
            if (!visible_block(coord + glm::ivec3(0, 0, 1))) {
                vert.packed_residue_face_bits |= pos_z_face_bit;
            }

            // Add this voxel only if it's visible.
            if ((vert.packed_residue_face_bits & all_face_bits) != 0) {
                assert(mesh.vert_count < chunk_max_verts);
                mesh.verts[mesh.vert_count++] = vert;
            }
        };

        mesh.vert_count = 0;

        for (int z = 0; z < chunk_size; ++z) {
            for (int y = 0; y < chunk_size; ++y) {
                for (int x = 0; x < chunk_size; ++x) {
                    visit_voxel(glm::ivec3(x, y, z));
                }
            }
        }
    }

    // Fill the given MeshEntry with mesh data for the given chunk
    // group from the given world.
    static void replace(PositionedChunkGroup& pcg,
                        VoxelWorld& world,
                        MeshEntry* entry)
    {
        entry->world_id = world.id();
        entry->group_coord = group_coord(pcg);

        if (entry->vbo_name == 0) {
            glGenBuffers(1, &entry->vbo_name);
            glBindBuffer(GL_ARRAY_BUFFER, entry->vbo_name);
            glBufferStorage(GL_ARRAY_BUFFER,
                            MeshEntry::vbo_bytes,
                            nullptr,
                            GL_DYNAMIC_STORAGE_BIT);
            PANIC_IF_GL_ERROR;
        }

        // Parasitically depend on the update, using the AlwaysDirty
        // parameter to get every chunk's mesh copied onto the new (or
        // reused) VBO even if not marked dirty.
        update<true>(pcg, world, entry, false);
    }

    // Update the given mesh entry with new data from dirty chunks,
    // except report failure (return false) if read_only is true
    // but dirty data needed to be re-uploaded anyway.
    //
    // Requires that the MeshEntry corresponds to the given chunk
    // group of the given world.
    template <bool AlwaysDirty=false>
    static bool update(PositionedChunkGroup& pcg,
                       VoxelWorld& world,
                       MeshEntry* entry,
                       bool read_only)
    {
        assert(entry->world_id == world.id());
        assert(entry->group_coord == group_coord(pcg));

        auto vbo_name = entry->vbo_name;
        assert(vbo_name != 0);

        // Recompute and reupload the mesh for any dirty chunk, and
        // reset the chunk dirty flag.
        for (int z = 0; z < edge_chunks; ++z) {
            for (int y = 0; y < edge_chunks; ++y) {
                for (int x = 0; x < edge_chunks; ++x) {
                    Chunk& chunk = group(pcg).chunk_array[z][y][x];
                    if (!AlwaysDirty and !chunk.mesh_dirty) continue;
                    if (chunk.total_visible == 0) continue;

                    if (read_only) return false;

                    ChunkMesh& mesh = entry->mesh_array[z][y][x];
                    glm::ivec3 residue = glm::ivec3(x,y,z) * chunk_size;
                    fill_mesh_verts(mesh, chunk, residue);
                    glNamedBufferSubData(vbo_name,
                                         entry->byte_offset(x, y, z),
                                         sizeof mesh.verts[0] * mesh.vert_count,
                                         mesh.verts);
                    chunk.mesh_dirty = false;
                }
            }
        }
        PANIC_IF_GL_ERROR;
        return true;
    }

    // Render, to the current framebuffer, chunks near the camera
    // using the conventional mesh-based algorithm.
    //
    // Each MeshEntry contains a list of visible voxels: their residue
    // coordinates, colors, and which of their 6 faces are visible.
    // My plan is to use instanced rendering to draw a bunch of unit
    // cubes in the correct locations, using degenerate triangles to
    // hide the hidden faces.
    void render_world_mesh_step() noexcept
    {
        MeshStore& store = camera.get_mesh_store();
        store.stats_begin_frame();

        glm::mat4 residue_vp_matrix = camera.get_residue_vp();
        glm::vec3 eye_residue;
        glm::ivec3 eye_group;
        camera.get_eye(&eye_group, &eye_residue);

        static GLuint vao = 0;
        static GLuint program_id;
        static GLint mvp_matrix_idx;

        // Position of camera eye relative to the origin of the group
        // (origin == group_size times the group coordinate).
        static GLint eye_relative_group_origin_id;

        static GLint far_plane_squared_id;
        static GLint fog_enabled_id;
        static GLint black_fog_id;

        if (vao == 0) {
            program_id = make_program({
                "mesh.vert", "mesh.frag", "fog_border.frag" });
            mvp_matrix_idx = glGetUniformLocation(program_id, "mvp_matrix");
            assert(mvp_matrix_idx >= 0);
            eye_relative_group_origin_id = glGetUniformLocation(program_id,
                "eye_relative_group_origin");
            assert(eye_relative_group_origin_id >= 0);
            far_plane_squared_id = glGetUniformLocation(program_id,
                "far_plane_squared");
            assert(far_plane_squared_id >= 0);
            fog_enabled_id = glGetUniformLocation(program_id,
                "fog_enabled");
            assert(fog_enabled_id >= 0);
            black_fog_id = glGetUniformLocation(program_id, "black_fog");
            assert(black_fog_id >= 0);

            glGenVertexArrays(1, &vao);
            glBindVertexArray(vao);
            PANIC_IF_GL_ERROR;

            // The above was all boilerplate crap; this is where I
            // bind once and forever the unit box buffers to the VAO.
            GLuint buffers[2];
            glGenBuffers(2, buffers);
            glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[1]);

            glBufferStorage(GL_ARRAY_BUFFER,
                            sizeof voxel_unit_box_vertices,
                            voxel_unit_box_vertices,
                            0);
            glBufferStorage(GL_ELEMENT_ARRAY_BUFFER,
                            sizeof voxel_unit_box_elements,
                            voxel_unit_box_elements,
                            0);

            glVertexAttribPointer(
                unit_box_vertex_idx,
                3,
                GL_FLOAT,
                false,
                sizeof(VoxelUnitBoxVertex),
                (void*) offsetof(VoxelUnitBoxVertex, x));
            glEnableVertexAttribArray(unit_box_vertex_idx);

            glVertexAttribIPointer(
                unit_box_face_bit_idx,
                1,
                GL_UNSIGNED_INT,
                sizeof(VoxelUnitBoxVertex),
                (void*) offsetof(VoxelUnitBoxVertex, face_bit));
            glEnableVertexAttribArray(unit_box_face_bit_idx);

            glVertexAttribPointer(
                unit_box_uv_idx,
                2,
                GL_FLOAT,
                false,
                sizeof(VoxelUnitBoxVertex),
                (void*) offsetof(VoxelUnitBoxVertex, u));
            glEnableVertexAttribArray(unit_box_uv_idx);

            PANIC_IF_GL_ERROR;
        }

        // If we're not using D-tier Intel drivers, this should
        // restore the element array binding above...
        glBindVertexArray(vao);

        glUseProgram(program_id);
        auto far_plane = camera.get_far_plane();
        glUniform1i(far_plane_squared_id, far_plane * far_plane);
        glUniform1i(fog_enabled_id, camera.get_fog());
        glUniform1i(black_fog_id, camera.use_black_fog());
        PANIC_IF_GL_ERROR;
        StoreRequest request;
        unsigned drawn_group_count = 0;
        unsigned drawn_chunk_count = 0;

        auto draw_group = [&] (PositionedChunkGroup& pcg)
        {
            if (group(pcg).total_visible == 0) return;
            MeshEntry* entry = store.request<Renderer>(pcg, world, &request);

            assert(entry->vbo_name != 0);
            glBindBuffer(GL_ARRAY_BUFFER, entry->vbo_name);

            // Bind instanced (per-voxel) attributes.
            // TODO: Maybe make one VAO per MeshEntry?
            glVertexAttribIPointer(
                packed_vertex_idx,
                1,
                GL_UNSIGNED_INT,
                sizeof(MeshVoxelVertex),
                (void*) offsetof(MeshVoxelVertex, packed_residue_face_bits));
            glEnableVertexAttribArray(packed_vertex_idx);
            glVertexAttribDivisor(packed_vertex_idx, 1);

            glVertexAttribIPointer(
                packed_color_idx,
                1,
                GL_UNSIGNED_INT,
                sizeof(MeshVoxelVertex),
                (void*) offsetof(MeshVoxelVertex, packed_color));
            glEnableVertexAttribArray(packed_color_idx);
            glVertexAttribDivisor(packed_color_idx, 1);
            PANIC_IF_GL_ERROR;

            // The view matrix only takes into account the eye's
            // residue coordinate, so the model position of the group
            // actually needs to be shifted by the eye's group coord.
            glm::vec3 model_offset = glm::vec3(group_coord(pcg) - eye_group)
                                   * float(group_size);
            glm::mat4 m = glm::translate(glm::mat4(1.0f), model_offset);
            glm::mat4 mvp = residue_vp_matrix * m;
            glUniformMatrix4fv(mvp_matrix_idx, 1, 0, &mvp[0][0]);

            // Similarly, the eye residue needs to be shifted by the
            // group's position.
            glm::vec3 eye_relative_group_origin = eye_residue - model_offset;
            glUniform3fv(eye_relative_group_origin_id, 1,
                &eye_relative_group_origin[0]);

            auto gc = group_coord(pcg);
            for (int z = 0; z < edge_chunks; ++z) {
                for (int y = 0; y < edge_chunks; ++y) {
                    for (int x = 0; x < edge_chunks; ++x) {
                        Chunk& chunk = group(pcg).chunk_array[z][y][x];
                        if (decide_chunk(gc, chunk) != draw_mesh) continue;

                        const ChunkMesh& mesh = entry->mesh_array[z][y][x];
                        if (mesh.vert_count == 0) continue;

                        // Longest function name EVER. I need it because
                        // every chunk's data is at a different offset
                        // within the VBO for instanced data (position/color).
                        glDrawElementsInstancedBaseVertexBaseInstance(
                            GL_TRIANGLES,
                            36, // 12 triangles, 3 verts each.
                            GL_UNSIGNED_SHORT,
                            nullptr,
                            mesh.vert_count,
                            0, // base vertex is always 0, for the unit box.
                            entry->vert_offset(x, y, z)
                            // ^^^ Instance offset, depends on chunk.
                        );
                        ++drawn_chunk_count;
                    }
                }
            }
            ++drawn_group_count;
            PANIC_IF_GL_ERROR;
        };

        float squared_thresh = camera.get_raycast_threshold();
        squared_thresh *= squared_thresh;

        // TODO: Avoid visiting too-far-away chunk groups.
        for (PositionedChunkGroup& pcg : world.group_map) {
            float min_squared_dist;
            bool cull = cull_group(pcg, &min_squared_dist);
            // Skip culled chunk groups or those so far away that they
            // cannot possibly contain chunks near enough to be drawn
            // using the mesh renderer.
            if (cull or min_squared_dist >= squared_thresh) continue;
            draw_group(pcg);
        }

        glBindVertexArray(0);

        if (evict_stats_debug) {
            fprintf(stderr, "\x1b[36mMeshStore stats:\x1b[0m\n");
            store.print_stats_end_frame();
            fprintf(stderr, "Mesh-drawn chunks: %u    groups: %u\n",
                drawn_chunk_count, drawn_group_count);
        }
    }

    // Now time to write the AABB-raycasting renderer. Here goes...

    // Basic idea for uploading voxel data to the GPU; maybe try to
    // explain better later if it works. Chunk groups are represented
    // by 3D textures stored in a LRU-eviction set-associative cache
    // (RaycastStore). When the texture is in the RaycastStore, its
    // contents are considered unchangable.
    //
    // If a new chunk group needs to be written into the RaycastStore,
    // or is dirty and needs to be updated, its contents are first
    // written to a StagingBuffer. This consists of a
    // persistent-mapped SSBO and a texture that is the target of a
    // compute shader using image load-store to transfer data from the
    // SSBO to the texture.
    //
    // When the compute shader is done, the filled texture is swapped
    // into the right place in RaycastStore (the swapped-out texture
    // is now usable for staging).
    //
    // This sounded like madness but Christoph Kubisch of Nvidia
    // convinced me to give it a try (the tiling of 3D textures is
    // really performant compared to a home-made SSBO).

    // Data flow of new or updated chunk group being written to a 3D texture:
    //
    // 1. Add chunk group to the queue of chunk groups waiting to be written.
    //
    // 2. When at the head of the queue (or near enough), get a staging
    //    buffer from the write_staging_buffers list. Copy into the
    //    persistent-mapped staging buffer. Note: memory barrier at (5) makes
    //    this write safe.
    //
    //    This is also when we mark the chunks as non-dirty, as at this
    //    point we're committed to updating the GPU storage.
    //
    // 3. Flush SSBO and dispatch compute.
    //
    // 4. Swap the write_staging_buffers and read_staging_buffers. Wait a frame.
    //
    // 5. Issue a memory barrier. At this points the read_staging_buffers
    //    have textures that are globally visible.
    //
    // 6. Swap the textures from read_staging_buffers into the correct places
    //    in the RaycastStore, and draw the modified chunk groups. At this
    //    point the texture that used to be in the RaycastStore (if any)
    //    is now used for staging (i.e. it's allowed to be modified).

    // Pseudocode:
    //
    // Mark all write staging buffers as non used.
    //
    // For all in queue (of world_id and group_coord pairs):
    //   Reserve unused write staging buffer if possible.
    //   Lookup chunk group in world and write to staging SSBO. (2)
    //   At this point also upload chunk AABBs
    //     (so they always match actually-uploaded data)
    //
    // For all chunk groups in the world visible in the view frustum:
    //   If dirty/not found in RaycastStore (new chunk group)
    //     Reserve write staging buffer if possible and write to SSBO (2)
    //       Also upload chunk AABBs as before.
    //     Otherwise add the world_id and group_coord to the queue (1)
    //
    //   If there is a matching staging buffer in read_staging_buffers
    //     Swap into RaycastStore.
    //     Add to list of chunk groups to draw after memory barrier. (6)
    //   Otherwise, add to before memory barrier list.
    //
    //   Note: For efficency, above pseudocode actually is done by
    //   scattering stuff in read_staging_buffers into the correct
    //   location in RaycastStore and marking a "memory barrier
    //   needed" bit.
    //
    //   Draw chunk groups in before memory barrier list.
    //   Issue memory barrier. (5)
    //   Draw chunk groups in after memory barrier list.
    //
    //   For all in-use StagingBuffers in write_staging_buffers:
    //     Dispatch compute shader. (3)
    //
    //   Swap write_staging_buffers and read_staging_buffers. (4)


    // Functions required by BaseStore. Since I'm revamping voxel upload
    // using the SSBO-approach above, I just provide dummy functions for
    // the BaseStore. So this design didn't age well...
    static void replace(PositionedChunkGroup& pcg,
                        VoxelWorld& world,
                        RaycastEntry* entry)
    {
        entry->world_id = world.id();
        entry->group_coord = group_coord(pcg);
        entry->should_draw = false;
    }

    // Another dummy function.
    static bool update(PositionedChunkGroup&,
                       VoxelWorld&,
                       RaycastEntry*,
                       bool)
    {
        return true;
    }

  private:
    // Mark all write staging buffers as unused.
    void mark_write_staging_buffers_unused()
    {
        RaycastStore& store = camera.get_raycast_store();
        for (StagingBuffer& sb : store.write_staging_buffers) {
            sb.world_id = 0;
        }
    }

    // Handle step 2 of the data flow: try to get a staging buffer
    // (return true iff this is done) and fill it with data from the
    // given PositionedChunkGroup (implicitly from the world being
    // rendered).
    //
    // We can clear the dirty flags exactly when this is done
    // successfully.
    bool try_upload_to_staging(PositionedChunkGroup& pcg)
    {
        RaycastStore& store = camera.get_raycast_store();
        StagingBuffer* p_sb = nullptr;
        for (StagingBuffer& sb : store.write_staging_buffers) {
            if (sb.world_id == 0) {
                p_sb = &sb;
                break;
            }
        }

        if (p_sb == nullptr) return false;

        // Label this StagingBuffer as ours.
        StagingBuffer& sb = *p_sb;
        sb.world_id = world.id();
        sb.group_coord = group_coord(pcg);
        sb.make_ready();

        // Fill in the texels and AABBs.
        ChunkGroup& cg = group(pcg);

        for (int z = 0; z < group_size; ++z) {
            auto& chunk_plane = cg.chunk_array[z >> chunk_shift];
            auto zl = z & (chunk_size - 1);

            for (int y = 0; y < group_size; ++y) {
                auto& chunk_row = chunk_plane[y >> chunk_shift];
                auto yl = y & (chunk_size - 1);

                for (int x = 0; x < group_size; ++x) {
                    Chunk& chunk = chunk_row[x >> chunk_shift];
                    auto xl = x & (chunk_size - 1);
                    chunk.texture_dirty = false;
                    Voxel voxel = chunk.voxel_array[zl][yl][xl];
                    sb.mapped_ssbo->voxel_colors[z][y][x] =
                        to_packed_color(voxel);
                }
            }
        }
        cg.dirty = false;

        // Fill in the AABBs.
        StagingBuffer::AABBs a;
        for (int z = 0; z < edge_chunks; ++z) {
            for (int y = 0; y < edge_chunks; ++y) {
                for (int x = 0; x < edge_chunks; ++x) {
                    Chunk& chunk = group(pcg).chunk_array[z][y][x];
                    glm::ivec3 low, high;
                    chunk.get_aabb(&low, &high);
                    a.aabb_array[z][y][x] = PackedAABB(low, high);
                }
            }
        }
        assert(sb.vbo_name != 0);
        glNamedBufferSubData(
            sb.vbo_name, 0, sizeof a.aabb_array, &a.aabb_array);
        PANIC_IF_GL_ERROR;

        return true;
    }

    // Pop stuff from the head (front) of the queue and start
    // uploading to the write staging buffers.
    void queue_to_staging_buffers()
    {
        RaycastStore& store = camera.get_raycast_store();

        auto it = store.queue.begin();

        for ( ; it != store.queue.end(); ++it) {
            auto world_id = it->world_id;
            auto group_coord = it->group_coord;

            if (world_id != world.id()) continue;

            auto pcg_it = world.group_map.find(group_coord);
            if (pcg_it == world.group_map.end()) continue;
            PositionedChunkGroup& pcg = *pcg_it;

            bool success = try_upload_to_staging(pcg);
            if (!success) break;
        }

        // Remove from queue up until the last successfully uploaded
        // (or discarded) chunk group.
        store.queue.erase(store.queue.begin(), it);
    }

    // Take stuff in the read staging buffers and put it where it
    // needs to go in the RaycastStore, marking the new entries as
    // needing a memory barrier.
    void swap_read_staging_buffers_into_RaycastStore()
    {
        RaycastStore& store = camera.get_raycast_store();
        for (StagingBuffer& sb : store.read_staging_buffers) {
            if (sb.world_id == 0) continue;

            StoreRequest request;
            RaycastEntry* entry =
                store.request_no_update(sb.group_coord, world.id(), &request);

            assert(entry != nullptr);
            using std::swap;
            assert(entry->world_id == sb.world_id);
            sb.world_id = 0;
            assert(entry->group_coord == sb.group_coord);
            swap(entry->vbo_name, sb.vbo_name);
            swap(entry->texture_name, sb.texture_name);
            entry->needs_memory_barrier = true;
            entry->should_draw = true;
        }
    }

    // Dispatch the compute shader that transfers pixels from staging
    // buffer SSBOs to the actual texture.
    void dispatch_write_staging_buffers_compute_shaders()
    {
        static GLuint program_id = 0;

        if (program_id == 0) {
            program_id = make_program( { "staging.comp" } );
        }
        glUseProgram(program_id);

        // Hook up the SSBO binding to the SSBO index that will
        // be used to deliver voxel data.
        glShaderStorageBlockBinding(program_id,
            chunk_group_voxels_program_index,
            chunk_group_voxels_binding_index);
        // Hook up the out_image to the correct image unit.
        glUniform1i(staging_image_program_index, staging_image_unit);

        PANIC_IF_GL_ERROR;
        RaycastStore& store = camera.get_raycast_store();
        for (StagingBuffer& sb : store.write_staging_buffers) {
            if (sb.world_id == 0) continue;
            auto sz = group_size * group_size * group_size * sizeof(uint32_t);
            assert(sb.ssbo_name != 0);
            glFlushMappedNamedBufferRange(sb.ssbo_name, 0, sz);
            PANIC_IF_GL_ERROR;

            // Now that the buffer is flushed, we can bind the source buffer
            // and output image and dispatch the compute shader.
            glBindBufferBase(
                GL_SHADER_STORAGE_BUFFER,
                chunk_group_voxels_binding_index,
                sb.ssbo_name);
            glBindImageTexture(
                staging_image_unit,
                sb.texture_name,
                0, false, 0, GL_WRITE_ONLY, GL_RGBA8UI);
            glDispatchCompute(8, 1, 1);
            PANIC_IF_GL_ERROR;
        }
        // RaycastStore& store = camera.get_raycast_store();
        // for (StagingBuffer& sb : store.write_staging_buffers) {
        //     PANIC_IF_GL_ERROR;
        //     if (sb.world_id == 0) continue;
        //     auto sz = group_size * group_size * group_size * sizeof(uint32_t);
        //     assert(sb.ssbo_name != 0);
        //     glFlushMappedNamedBufferRange(sb.ssbo_name, 0, sz);
        //     PANIC_IF_GL_ERROR;
        //     glBindBuffer(GL_PIXEL_UNPACK_BUFFER, sb.ssbo_name);
        //     PANIC_IF_GL_ERROR;
        //     glTextureSubImage3D(
        //         sb.texture_name, 0, 0, 0, 0,
        //         group_size, group_size, group_size,
        //         GL_RGBA, GL_UNSIGNED_INT_8_8_8_8, nullptr);
        //     glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
        // }
        // PANIC_IF_GL_ERROR;
    }

  public:
    // Render, to the current framebuffer, chunks around the camera
    // using the AABB-raycast algorithm. Chunks that are near
    // enough to have been drawn using the mesh algorithm will
    // not be re-drawn.
    void render_world_raycast_step() noexcept
    {
        RaycastStore& store = camera.get_raycast_store();
        store.stats_begin_frame();

        // Resize the cache to a size suitable for the given render distance.
        auto far_plane = camera.get_far_plane();
        uint32_t modulus = uint32_t(std::max(4.0, ceil(far_plane / 128.0)));
        store.set_modulus(modulus, world.id());
        store.resize_victim_cache(0);

        // Count and limit the number of new chunk groups added to the
        // RaycastStore this frame.
        int remaining_new_chunk_groups_allowed =
            camera.get_max_frame_new_chunk_groups();

        // Collect all chunk groups needing to be raycast, and sort from
        // nearest to furthest (reverse painters). This fascilitates
        // the early depth test optimization.
        float squared_thresh = camera.get_far_plane();
        squared_thresh *= squared_thresh;
        std::vector<std::pair<float, PositionedChunkGroup*>> pcg_by_depth;

        // TODO: Avoid visiting far away chunk groups.
        unsigned drawn_group_count = 0;
        for (PositionedChunkGroup& pcg : world.group_map) {
            float min_squared_dist;
            bool cull = cull_group(pcg, &min_squared_dist);
            if (cull or min_squared_dist >= squared_thresh) continue;
            pcg_by_depth.emplace_back(min_squared_dist, &pcg);
            ++drawn_group_count;
        }

        if (!disable_zcull_sort) {
            std::sort(pcg_by_depth.begin(), pcg_by_depth.end());
        }
        std::vector<RaycastEntry*> draw_in_first_pass;
        std::vector<RaycastEntry*> draw_in_second_pass;

        // Implement the earlier pseudocode here.
        store.write_staging_buffers.resize(96);
        store.read_staging_buffers.resize(96);
        mark_write_staging_buffers_unused();
        queue_to_staging_buffers(); // (2)

        for (auto& pair : pcg_by_depth) {
            // Steps 1 and 2: prepare to add dirty chunk groups to staging.
            StoreRequest request;
            request.may_fail = remaining_new_chunk_groups_allowed <= 0;
            PositionedChunkGroup& pcg = *pair.second;
            RaycastEntry* entry = store.request<Renderer>(pcg, world, &request);

            if (entry == nullptr) continue;

            ChunkGroup& cg = group(pcg);
            if (cg.dirty or !request.was_in_store) {
                bool success = try_upload_to_staging(pcg); // (2)
                if (!success) {
                    // (1)
                    store.queue.push_back( { world.id(), group_coord(pcg) } );
                }
            }
            remaining_new_chunk_groups_allowed -= !request.was_in_store;
        }

        swap_read_staging_buffers_into_RaycastStore(); // (6)
        for (auto& pair : pcg_by_depth) {
            StoreRequest request;
            request.may_fail = true;
            PositionedChunkGroup& pcg = *pair.second;
            RaycastEntry* entry = store.request<Renderer>(pcg, world, &request);

            if (entry == nullptr) continue;
            if (!entry->should_draw) continue;

            if (entry->needs_memory_barrier) {
                draw_in_second_pass.push_back(entry);
            }
            else {
                draw_in_first_pass.push_back(entry);
            }
        }

        draw_raycast_entries(draw_in_first_pass, false);
        glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT
                      | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT); // (5)
        draw_raycast_entries(draw_in_second_pass, true); // (6)

        dispatch_write_staging_buffers_compute_shaders(); // (3)

        swap(store.read_staging_buffers, store.write_staging_buffers);

        if (evict_stats_debug) {
            fprintf(stderr, "\x1b[36mRaycastStore stats:\x1b[0m\n");
            store.print_stats_end_frame();
            fprintf(stderr, "Raycast groups: %u\n", drawn_group_count);
        }
    }

  private:
    // Draw the given list of RaycastEntry (they have the needed group
    // coord).  The should_have_memory_barrier_bit is for checking
    // purposes: it should be set to false before the memory barrier
    // was issued (and we check none of the RaycastEntries wanted a
    // barrier), and should be set to true after issue (in which case
    // we clear the bit to false as the barrier has been issued
    // already).
    void draw_raycast_entries(
        const std::vector<RaycastEntry*>& entries,
        bool should_have_memory_barrier_bit)
    {
        glm::mat4 residue_vp_matrix = camera.get_residue_vp();
        glm::vec3 eye_residue;
        glm::ivec3 eye_group;
        camera.get_eye(&eye_group, &eye_residue);
        auto far_plane = camera.get_far_plane();

        static GLuint vao = 0;
        static GLuint program_id;
        static GLuint vertex_buffer_id;
        static GLuint element_buffer_id;
        static GLint mvp_matrix_id;
        // Position of camera eye relative to the origin of the group
        // (origin == group_size times the group coordinate).
        static GLint eye_relative_group_origin_id;
        static GLint far_plane_squared_id;
        static GLint raycast_thresh_squared_id;
        static GLint fog_enabled_id;
        static GLint black_fog_id;
        static GLint chunk_debug_id;
        static GLint chunk_group_texture_id;

        if (vao == 0) {
            PANIC_IF_GL_ERROR;
            program_id = make_program({
                "raycast.vert", "raycast.frag",
                "fog_border.frag", "read_group_voxel.frag" });
            mvp_matrix_id = glGetUniformLocation(program_id, "mvp_matrix");
            assert(mvp_matrix_id >= 0);
            eye_relative_group_origin_id = glGetUniformLocation(program_id,
                "eye_relative_group_origin");
            assert(eye_relative_group_origin_id >= 0);
            chunk_debug_id = glGetUniformLocation(program_id, "chunk_debug");
            assert(chunk_debug_id >= 0);
            far_plane_squared_id = glGetUniformLocation(program_id,
                "far_plane_squared");
            assert(far_plane_squared_id >= 0);
            raycast_thresh_squared_id = glGetUniformLocation(program_id,
                "raycast_thresh_squared");
            assert(raycast_thresh_squared_id >= 0);
            fog_enabled_id = glGetUniformLocation(program_id,
                "fog_enabled");
            assert(fog_enabled_id >= 0);
            black_fog_id = glGetUniformLocation(program_id, "black_fog");
            assert(black_fog_id >= 0);
            chunk_group_texture_id =
                glGetUniformLocation(program_id, "chunk_group_texture");
            assert(chunk_group_texture_id >= 0);

            glGenVertexArrays(1, &vao);
            glBindVertexArray(vao);

            // The binding points for the non-instanced unit cube
            // vertices, and the element buffer binding, will be
            // stored forever in the VAO.
            glGenBuffers(1, &vertex_buffer_id);
            glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_id);
            glBufferStorage(
                GL_ARRAY_BUFFER, sizeof unit_box_vertices,
                unit_box_vertices, 0);

            glGenBuffers(1, &element_buffer_id);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer_id);
            glBufferStorage(
                GL_ELEMENT_ARRAY_BUFFER, sizeof unit_box_elements,
                unit_box_elements, 0);

            glVertexAttribPointer(
                unit_box_vertex_idx,
                3,
                GL_FLOAT,
                false,
                sizeof(float) * 6,
                (void*)0);
            glEnableVertexAttribArray(unit_box_vertex_idx);

            glVertexAttribPointer(
                unit_box_normal_idx,
                3,
                GL_FLOAT,
                false,
                sizeof(float) * 6,
                (void*)12);
            glEnableVertexAttribArray(unit_box_normal_idx);
        }

        glBindVertexArray(vao);
        glUseProgram(program_id);

        // Set uniforms unchanged per-chunk-group.
        glUniform1i(chunk_debug_id, chunk_debug);
        glUniform1i(fog_enabled_id, camera.get_fog());
        glUniform1i(black_fog_id, camera.use_black_fog());
        glUniform1i(far_plane_squared_id, far_plane * far_plane);
        auto raycast_thr = camera.get_raycast_threshold();
        glUniform1i(raycast_thresh_squared_id, raycast_thr * raycast_thr);

        // I'll use texture 0 for the chunk group texture.
        glActiveTexture(GL_TEXTURE0);
        glUniform1i(chunk_group_texture_id, 0);

        PANIC_IF_GL_ERROR;

        // My plan is to use instanced rendering to draw the chunk AABBs
        // of this chunk group. The base box is a 1x1x1 unit box, which
        // is stretched and repositioned in the vertex shader to the
        // true AABB.
        for (RaycastEntry* entry : entries) {
            assert(entry->texture_name != 0);
            assert(entry->vbo_name != 0);
            if (entry->needs_memory_barrier != should_have_memory_barrier_bit) {
                fprintf(stderr, "May have missed memory barrier...\n");
            }
            entry->needs_memory_barrier = false;

            // The view matrix only takes into account the eye's
            // residue coordinate, so the model position of the group
            // actually needs to be shifted by the eye's group coord.
            glm::vec3 model_offset = glm::vec3(entry->group_coord - eye_group)
                                   * float(group_size);
            glm::mat4 m = glm::translate(glm::mat4(1.0f), model_offset);
            glm::mat4 mvp = residue_vp_matrix * m;
            glUniformMatrix4fv(mvp_matrix_id, 1, 0, &mvp[0][0]);
            PANIC_IF_GL_ERROR;

            // Similarly, the eye residue needs to be shifted by the
            // group's position.
            glm::vec3 eye_relative_group_origin = eye_residue - model_offset;
            glUniform3fv(eye_relative_group_origin_id, 1,
                &eye_relative_group_origin[0]);

            // The unit box vertex attribs should already be bound by VAO.
            //
            // Get the instanced vertex attribs going (i.e. bind the
            // AABB residue coords, packed as integers). I should
            // probably make a VAO per chunk group in the
            // RaycastEntry.
            glBindBuffer(GL_ARRAY_BUFFER, entry->vbo_name);
            glVertexAttribIPointer(
                packed_aabb_low_idx,
                1,
                GL_UNSIGNED_INT,
                sizeof(PackedAABB),
                (void*) offsetof(PackedAABB, packed_low));
            glEnableVertexAttribArray(packed_aabb_low_idx);
            glVertexAttribDivisor(packed_aabb_low_idx, 1);

            glVertexAttribIPointer(
                packed_aabb_high_idx,
                1,
                GL_UNSIGNED_INT,
                sizeof(PackedAABB),
                (void*) offsetof(PackedAABB, packed_high));
            glEnableVertexAttribArray(packed_aabb_high_idx);
            glVertexAttribDivisor(packed_aabb_high_idx, 1);
            PANIC_IF_GL_ERROR;

            // Bind the relevant chunk group texture.
            glBindTexture(GL_TEXTURE_3D, entry->texture_name);

            // Draw all edge_chunks^3 chunks in the chunk group.
            auto instances = edge_chunks * edge_chunks * edge_chunks;
            glDrawElementsInstanced(
                GL_TRIANGLES, 36, GL_UNSIGNED_SHORT, nullptr, instances);
            PANIC_IF_GL_ERROR;
        }

        PANIC_IF_GL_ERROR;
        glBindVertexArray(0);
    }
};

void render_world_mesh_step(VoxelWorld& world, Camera& camera)
{
    Renderer(world, camera).render_world_mesh_step();
}

void render_world_raycast_step(VoxelWorld& world, Camera& camera)
{
    Renderer(world, camera).render_world_raycast_step();
}

// Render the background. This uses a shader hard-wired to draw a
// full-screen rectangle.
void render_background(Camera& camera)
{
    static GLuint vao = 0;
    static GLuint program_id;
    static GLint inverse_vp_id;
    static GLint eye_world_position_id;
    static GLint black_fog_id;

    if (vao == 0) {
        glGenVertexArrays(1, &vao);
        program_id = make_program(
            { "background.vert", "background.frag", "fog_border.frag" } );

        inverse_vp_id = glGetUniformLocation(program_id,
            "inverse_vp");
        assert(inverse_vp_id >= 0);
        eye_world_position_id = glGetUniformLocation(program_id,
            "eye_world_position");
        assert(eye_world_position_id >= 0);
        black_fog_id = glGetUniformLocation(program_id,
            "black_fog");
        assert(black_fog_id >= 0);
        PANIC_IF_GL_ERROR;
    }

    glm::vec3 eye_residue;
    camera.get_eye(nullptr, &eye_residue);
    glm::mat4 inverse_vp = glm::inverse(camera.get_residue_vp());

    glBindVertexArray(vao);
    glUseProgram(program_id);
    glUniformMatrix4fv(inverse_vp_id, 1, 0, &inverse_vp[0][0]);
    glUniform3fv(eye_world_position_id, 1, &eye_residue[0]);
    glUniform1i(black_fog_id, camera.use_black_fog());
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    glBindVertexArray(0);
    PANIC_IF_GL_ERROR;
};

MeshStore* new_mesh_store()
{
    return new MeshStore;
}

void delete_mesh_store(MeshStore* mesh_store)
{
    delete mesh_store;
}

RaycastStore* new_raycast_store()
{
    return new RaycastStore;
}

void delete_raycast_store(RaycastStore* raycast_store)
{
    delete raycast_store;
}

void viewport(int x, int y)
{
    glViewport(0, 0, x, y);
    PANIC_IF_GL_ERROR;
};

void gl_first_time_setup()
{
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
}

void gl_clear()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
}

static GLuint f32_depth_framebuffer = 0;
static GLuint f32_depth_texture = 0;
static GLuint f32_depth_renderbuffer = 0;
static int f32_depth_framebuffer_x, f32_depth_framebuffer_y = 0;

// Used to detect requested screen size changes.
static int prev_screen_x = 0, prev_screen_y = 0, prev_target_fragments = 0;

// Ratio of screen size (length / width) to framebuffer size.
static int framebuffer_divisor = 1;


void bind_global_f32_depth_framebuffer(
    int screen_x, int screen_y, int target_fragments)
{
    int want_framebuffer_x = f32_depth_framebuffer_x;
    int want_framebuffer_y = f32_depth_framebuffer_y;

    if (target_fragments <= 0) {
        want_framebuffer_x = screen_x;
        want_framebuffer_y = screen_y;
        framebuffer_divisor = 1;
    }

    // Recalculate the actual framebuffer size if needed.
    else if (screen_x != prev_screen_x
         or screen_y != prev_screen_y
         or prev_target_fragments != target_fragments) {

        for (framebuffer_divisor = 1;
             framebuffer_divisor <= 16;
             ++framebuffer_divisor) {

            want_framebuffer_x = screen_x / framebuffer_divisor;
            want_framebuffer_y = screen_y / framebuffer_divisor;
            if (want_framebuffer_x * want_framebuffer_y <= target_fragments) {
                break;
            }
        }
    }

    prev_screen_x = screen_x;
    prev_screen_y = screen_y;
    prev_target_fragments = target_fragments;

    // Destroy the framebuffer and its attachments if the screen size changed.
    if (want_framebuffer_x != f32_depth_framebuffer_x
    or want_framebuffer_y != f32_depth_framebuffer_y) {
        if (f32_depth_framebuffer != 0) {
            fprintf(stderr, "Destroying old %i x %i framebuffer.\n",
                f32_depth_framebuffer_x, f32_depth_framebuffer_y);
            glDeleteFramebuffers(1, &f32_depth_framebuffer);
            glDeleteTextures(1, &f32_depth_texture);
            glDeleteRenderbuffers(1, &f32_depth_renderbuffer);
        }

        f32_depth_framebuffer = 0;
        f32_depth_framebuffer_x = want_framebuffer_x;
        f32_depth_framebuffer_y = want_framebuffer_y;
        PANIC_IF_GL_ERROR;
    }

    // (Re-) create the framebuffer if needed, and bind it.
    if (f32_depth_framebuffer == 0) {
        // Create framebuffer.
        fprintf(stderr, "Creating %i x %i framebuffer.\n",
                f32_depth_framebuffer_x, f32_depth_framebuffer_y);
        fprintf(stderr, "Downsampling by %i x %i.\n",
                framebuffer_divisor, framebuffer_divisor);

        glCreateFramebuffers(1, &f32_depth_framebuffer);
        glBindFramebuffer(GL_FRAMEBUFFER, f32_depth_framebuffer);
        PANIC_IF_GL_ERROR;

        // Add depth buffer.
        glCreateRenderbuffers(1, &f32_depth_renderbuffer);
        glNamedRenderbufferStorage(f32_depth_renderbuffer,
                                   GL_DEPTH_COMPONENT32F,
                                   f32_depth_framebuffer_x,
                                   f32_depth_framebuffer_y);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER,
                                  GL_DEPTH_ATTACHMENT,
                                  GL_RENDERBUFFER,
                                  f32_depth_renderbuffer);
        PANIC_IF_GL_ERROR;

        // Add color buffer.
        glCreateTextures(GL_TEXTURE_2D, 1, &f32_depth_texture);
        glTextureParameteri(f32_depth_texture, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTextureParameteri(f32_depth_texture, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTextureParameteri(f32_depth_texture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTextureParameteri(f32_depth_texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTextureStorage2D(f32_depth_texture, 1, GL_RGBA8,
                           f32_depth_framebuffer_x, f32_depth_framebuffer_y);
        glFramebufferTexture2D(GL_FRAMEBUFFER,
                               GL_COLOR_ATTACHMENT0,
                               GL_TEXTURE_2D,
                               f32_depth_texture,
                               0);
        PANIC_IF_GL_ERROR;

        // Wire up the only color output.
        GLenum tmp = GL_COLOR_ATTACHMENT0;
        glDrawBuffers(1, &tmp);

        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        {
            panic("Incomplete 32-bit depth framebuffer: "
                 + std::to_string(glCheckFramebufferStatus(GL_FRAMEBUFFER)));
        }
    }
    else {
        glBindFramebuffer(GL_FRAMEBUFFER, f32_depth_framebuffer);
    }
    glViewport(0, 0, f32_depth_framebuffer_x, f32_depth_framebuffer_y);
    PANIC_IF_GL_ERROR;
}

void finish_global_f32_depth_framebuffer(int screen_x, int screen_y)
{
    assert(f32_depth_framebuffer != 0);
    glBlitNamedFramebuffer(
        f32_depth_framebuffer,
        0, // Write to window framebuffer
        0, 0,
        f32_depth_framebuffer_x, f32_depth_framebuffer_y,
        0, 0,
        screen_x, screen_y,
        GL_COLOR_BUFFER_BIT,
        GL_NEAREST);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    PANIC_IF_GL_ERROR;
}

} // end namespace
