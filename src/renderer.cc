// Implementation of the hybrid voxel renderer. I don't see how I can
// make much of this exception safe, so I wrap everything in a
// noexcept at the end. As usual with OpenGL, none of this is
// thread-safe either (but this could change later now that I'm
// factoring everything into threads and may get rid of all the darn
// static variables).

#include "myricube.hh"

#include <algorithm>
#include <atomic>
#include <memory>
#include <stdio.h>
#include <typeinfo>
#include <thread>
#include <utility>
#include <vector>

#include "camera.hh"
#include "chunk.hh"
#include "glad/glad.h"
#include "renderer.hh"
#include "shaders.hh"
#include "window.hh"

namespace myricube {

// "Temporary" I used to support the idea of changing the world being
// rendered but not anymore, and they were distinguished with id's.
// Since I stopped this, for now I just give the global world a bogus ID.
//
// Eventually I should rip this out, but world_id also has a side
// effect of identifying valid/invalid cache entries.
const uint64_t bogus_world_id = 1;

// Hacky debug variables.
bool chunk_debug = false;
bool evict_stats_debug = false;
bool disable_zcull_sort = false;

// (roughly) minimum distance from the camera that a chunk needs
// to be to switch from mesh to raycast graphics.
// Keep as int to avoid rounding errors in distance culling.
constexpr int raycast_threshold = 120;

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

    // Compute the AABB (in residue coordinates) of the given chunk.
    // The residue coordinates of the lower corner of the chunk must
    // also be given (as the AABB is relative to the chunk group's
    // origin, not the chunk itself).
    PackedAABB(const BinChunk& chunk, glm::ivec3 chunk_residue)
    {
        // {xyz}_array[n] has visible_bit set iff there's some visible
        // voxel in the chunk with n == the voxel's {xyz} coordinate.
        uint32_t x_array[chunk_size] = { 0 };
        uint32_t y_array[chunk_size] = { 0 };
        uint32_t z_array[chunk_size] = { 0 };

        uint32_t all_orrd = 0;

        for (size_t z = 0; z < chunk_size; ++z) {
            for (size_t y = 0; y < chunk_size; ++y) {
                for (size_t x = 0; x < chunk_size; ++x) {
                    uint32_t voxel = chunk.voxel_array[z][y][x];
                    all_orrd |= voxel;
                    x_array[x] |= voxel;
                    y_array[y] |= voxel;
                    z_array[z] |= voxel;
                }
            }
        }

        // Return trivial PackedAABB if chunk is empty.
        if ((all_orrd & visible_bit) == 0) {
            packed_low = 0;
            packed_high = 0;
            return;
        }

        // Otherwise, we actually have to calculate the AABB. Do it
        // in the coordinate system of the chunk first.
        glm::ivec3 chunk_aabb_low(chunk_size, chunk_size, chunk_size);
        glm::ivec3 chunk_aabb_high(0, 0, 0);

        for (int32_t i = 0; i < chunk_size; ++i) {
            if (x_array[i] & visible_bit) {
                chunk_aabb_low.x = glm::min(chunk_aabb_low.x, i);
                chunk_aabb_high.x = glm::max(chunk_aabb_low.x, i+1);
            }
            if (y_array[i] & visible_bit) {
                chunk_aabb_low.y = glm::min(chunk_aabb_low.y, i);
                chunk_aabb_high.y = glm::max(chunk_aabb_low.y, i+1);
            }
            if (z_array[i] & visible_bit) {
                chunk_aabb_low.z = glm::min(chunk_aabb_low.z, i);
                chunk_aabb_high.z = glm::max(chunk_aabb_low.z, i+1);
            }
        }

        // Reposition computed AABB by the chunk's residue coordinate.
        *this = PackedAABB(
            chunk_aabb_low + chunk_residue,
            chunk_aabb_high + chunk_residue);
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
// single chunk. The AABB is also needed for now (decide_chunk needs
// this info).
struct ChunkMesh
{
    MeshVoxelVertex verts[chunk_max_verts];
    size_t vert_count = 0;
    PackedAABB aabb;
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
    // If possible, return a valid, updated Entry for the chunk group
    // with the given group coordinate. See StoreRequest for control
    // parameters.
    //
    // This requires a template parameter providing static functions:
    //
    // replace(glm::ivec3 group_coord, Entry*)
    //
    //     Fill the Entry with data for this new chunk group,
    //     initializing OpenGL resources if needed.
    //
    // bool update(glm::ivec3 group_coord, Entry*, bool read_only)
    //
    //     Assuming replace was already called with identical
    //     arguments, above, just re-upload any changed data, unless
    //     read_only is true. Return true iff the Entry is
    //     up-to-date (i.e. return false only when read_only was
    //     true, but some changed data needed to be re-uploaded).
    template <typename EntryFiller>
    Entry* request(EntryFiller& filler,
                   glm::ivec3 group_coord,
                   StoreRequest* request)
    {
        assert(!request->read_only or request->may_fail);

        bool valid;
        Entry* entry =
            cached_location(&valid,
                            group_coord,
                            bogus_world_id,
                            request->may_fail);
        request->was_in_store = valid;

        // If the Entry in the array already corresponds to the given
        // chunk group; just update it if we're allowed to (deal with
        // dirty chunks) and return.
        if (valid) {
            bool success =
                filler.update(group_coord, entry, request->read_only);
            assert(success or request->read_only);
            return success ? entry : nullptr;
        }

        // Otherwise, either fail if allowed, or evict and replace a
        // cache entry with data for the new chunk group.
        if (request->may_fail) return nullptr;

        filler.replace(group_coord, entry);
        assert(entry->world_id == bogus_world_id);
        assert(entry->group_coord == group_coord);
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
    //
    // (Actually, I'm now planning to fundamentally re-think my device
    // management scheme to allow queueing work and having worker
    // threads process the queue, in preparation for learning Vulkan).
    Entry* request_no_update(glm::ivec3 group_coord,
                             StoreRequest* request)
    {
        assert(!request->read_only or request->may_fail);

        bool valid;
        Entry* entry =
            cached_location(&valid,
                            group_coord,
                            bogus_world_id,
                            request->may_fail);
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

        entry->world_id = bogus_world_id;
        entry->group_coord = group_coord;
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
    void set_modulus(uint32_t new_modulus) noexcept
    {
        uint64_t world_id = bogus_world_id;
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

// Renderer class. Stores state that is used by a voxel rendering
// thread running for the lifetime of the Renderer.
class Renderer
{
    std::shared_ptr<Window> window_ptr;
    std::shared_ptr<SyncCamera> sync_camera_ptr;
    WorldHandle world_handle;
    ViewWorldCache world_cache;
    const BinChunkGroup* current_chunk_group = nullptr;

    // Updated from SyncCamera per frame.
    CameraTransforms tr;

    // Device storage for mesh and raycast rendering respectively.
    MeshStore mesh_store;
    RaycastStore raycast_store;

    // Used for communicating with the render thread.
    std::atomic<bool> thread_exit_flag { false };
    std::atomic<double> fps { 0 };
    std::atomic<double> frame_time { 0 };

    // Declare thread LAST so that thread starts with Renderer fully initialized.
    std::thread thread;
  public:
    Renderer(
        std::shared_ptr<Window> window_,
        const WorldHandle& world_,
        std::shared_ptr<SyncCamera> sync_camera_) :

        window_ptr(std::move(window_)),
        sync_camera_ptr(std::move(sync_camera_)),
        world_handle(world_),
        world_cache(world_),
        thread(render_loop, this) { }

    ~Renderer()
    {
        thread_exit_flag.store(true);
        thread.join();
    }

    Renderer(Renderer&&) = delete;

    friend double get_fps(const Renderer& renderer);
    friend double get_frame_time(const Renderer& renderer);

  private:
    // Sets the current_chunk_group ptr to point to the named chunk group.
    // (nullptr if it doesn't exist). This seems oddly stateful for my
    // tastes but I did it anyway to avoid causing a dangling pointer.
    // (Subsequent world_cache.get_entry may overwrite ptr).
    //
    // I may re-think the design of this program at some point. It used to
    // be quite clean but the quality has eroded with time as I stapled-on
    // more optimizations. But for now I'm just having fun.
    void set_current_chunk_group(glm::ivec3 group_coord)
    {
        auto& entry = world_cache.get_entry(group_coord);
        current_chunk_group = entry.chunk_group_ptr.get();
    }

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
    bool cull_group(glm::ivec3 group_coord, float* squared_dist=nullptr)
    {
        glm::mat4 vp = tr.residue_vp_matrix;
        glm::ivec3 eye_group = tr.eye_group;
        glm::vec3 eye_residue = tr.eye_residue;

        // Never cull the group the eye is in (this avoids
        // pathological cases e.g. 7 corners behind eye and 1 in front
        // and out of the frustum).
        if (group_coord == eye_group) {
            if (squared_dist) *squared_dist = 0.0f;
            return false;
        }

        eye_residue = glm::floor(eye_residue); // to match decide_chunk.

        // Position of this chunk group relative to the group that the
        // eye is in, in voxel units.
        glm::vec3 low_corner(group_size * (group_coord - eye_group));

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
                     PackedAABB aabb)
    {
        glm::ivec3 aabb_low(aabb.low_x(), aabb.low_y(), aabb.low_z());
        glm::ivec3 aabb_high(aabb.high_x(), aabb.high_y(), aabb.high_z());

        glm::ivec3 eye_group = tr.eye_group;
        glm::vec3 eye_residue = tr.eye_residue;
        auto far_plane = tr.far_plane;
        auto raycast_thresh = raycast_threshold;

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
                                const BinChunk& chunk,
                                glm::ivec3 chunk_residue)
    {
        mesh.aabb = PackedAABB(chunk, chunk_residue);

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
            return chunk(coord) & visible_bit;
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

  public:
    // Fill the given MeshEntry with mesh data for the given chunk
    // group from the given world.
    void replace(glm::ivec3 group_coord, MeshEntry* entry)
    {
        entry->world_id = bogus_world_id;
        entry->group_coord = group_coord;

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
        update<true>(group_coord, entry, false);
    }

    // Update the given mesh entry with new data from dirty chunks in
    // the named chunk group, except report failure (return false) if
    // read_only is true but dirty data needed to be re-uploaded
    // anyway.
    //
    // Requires that the MeshEntry corresponds to the given chunk
    // group of the given world.
    template <bool AlwaysDirty=false>
    bool update(glm::ivec3 group_coord,
                       MeshEntry* entry,
                       bool read_only)
    {
        assert(entry->world_id == bogus_world_id);
        assert(entry->group_coord == group_coord);

        auto vbo_name = entry->vbo_name;
        assert(vbo_name != 0);

        // Load in the needed chunk group.
        set_current_chunk_group(group_coord);
        assert(current_chunk_group != nullptr);

        bool dirty = AlwaysDirty ||
            current_chunk_group->dirty_flags & renderer_mesh_dirty_flag;

        if (!dirty) return true;
        if (read_only) return false;

        // Recompute and reupload the mesh for every chunk and reset
        // the dirty flag BEFORE reading (handles data races better).
        current_chunk_group->dirty_flags &= ~renderer_mesh_dirty_flag;
        for (int z = 0; z < edge_chunks; ++z) {
            for (int y = 0; y < edge_chunks; ++y) {
                for (int x = 0; x < edge_chunks; ++x) {
                    const BinChunk& chunk = current_chunk_group
                                          ->chunk_array[z][y][x];
                    ChunkMesh& mesh = entry->mesh_array[z][y][x];
                    glm::ivec3 residue = glm::ivec3(x,y,z) * chunk_size;
                    fill_mesh_verts(mesh, chunk, residue);
                    glNamedBufferSubData(vbo_name,
                                         entry->byte_offset(x, y, z),
                                         sizeof mesh.verts[0] * mesh.vert_count,
                                         mesh.verts);
                }
            }
        }

        PANIC_IF_GL_ERROR;
        return true;
    }

  private:
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
        MeshStore& store = mesh_store;
        store.stats_begin_frame();

        glm::mat4 residue_vp_matrix = tr.residue_vp_matrix;
        glm::vec3 eye_residue = tr.eye_residue;
        glm::ivec3 eye_group = tr.eye_group;

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
        }

        glBindVertexArray(vao);

        glUseProgram(program_id);
        auto far_plane = tr.far_plane;
        glUniform1i(far_plane_squared_id, far_plane * far_plane);
        glUniform1i(fog_enabled_id, tr.use_fog);
        glUniform1i(black_fog_id, tr.use_black_fog);
        PANIC_IF_GL_ERROR;
        StoreRequest request;
        unsigned drawn_group_count = 0;
        unsigned drawn_chunk_count = 0;

        auto draw_group = [&] (glm::ivec3 group_coord)
        {
            MeshEntry* entry = store.request(*this, group_coord, &request);

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
            glm::vec3 model_offset = glm::vec3(group_coord - eye_group)
                                   * float(group_size);
            glm::mat4 m = glm::translate(glm::mat4(1.0f), model_offset);
            glm::mat4 mvp = residue_vp_matrix * m;
            glUniformMatrix4fv(mvp_matrix_idx, 1, 0, &mvp[0][0]);

            // Similarly, the eye residue needs to be shifted by the
            // group's position.
            glm::vec3 eye_relative_group_origin = eye_residue - model_offset;
            glUniform3fv(eye_relative_group_origin_id, 1,
                &eye_relative_group_origin[0]);

            for (int z = 0; z < edge_chunks; ++z) {
                for (int y = 0; y < edge_chunks; ++y) {
                    for (int x = 0; x < edge_chunks; ++x) {
                        PackedAABB aabb = entry->mesh_array[z][y][x].aabb;
                        if (decide_chunk(group_coord, aabb) != draw_mesh) {
                            continue;
                        }

                        const ChunkMesh& mesh = entry->mesh_array[z][y][x];
                        if (mesh.vert_count == 0) continue;

                        // Need to choose the base instance because
                        // every chunk's data is at a different offset
                        // within the VBO for instanced data (position/color).
                        //
                        // As requested by mesh.vert (includes
                        // built-in model of a cube), need to draw as
                        // GL_TRIANGLES with 36 vertices.
                        glDrawArraysInstancedBaseInstance(
                            GL_TRIANGLES,
                            0, 36,
                            mesh.vert_count,
                            entry->vert_offset(x, y, z));
                            // ^^^ Instance offset, depends on chunk.
                        ++drawn_chunk_count;
                    }
                }
            }
            ++drawn_group_count;
            PANIC_IF_GL_ERROR;
        };

        float squared_thresh = raycast_threshold * raycast_threshold;

        glm::ivec3 group_coord_low, group_coord_high;
        glm::dvec3 disp(tr.far_plane);
        split_coordinate(tr.eye - disp, &group_coord_low);
        split_coordinate(tr.eye + disp, &group_coord_high);

        for (int32_t zH = group_coord_low.z; zH <= group_coord_high.z; ++zH) {
        for (int32_t yH = group_coord_low.y; yH <= group_coord_high.y; ++yH) {
        for (int32_t xH = group_coord_low.x; xH <= group_coord_high.x; ++xH) {
            auto group_coord = glm::ivec3(xH, yH, zH);

            set_current_chunk_group(group_coord);
            if (current_chunk_group == nullptr) continue;

            float min_squared_dist;
            bool cull = cull_group(group_coord, &min_squared_dist);
            // Skip culled chunk groups or those so far away that they
            // cannot possibly contain chunks near enough to be drawn
            // using the mesh renderer.
            if (cull or min_squared_dist >= squared_thresh) continue;
            draw_group(group_coord);
        }
        }
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

  public:
    // Functions required by BaseStore. Since I'm revamping voxel upload
    // using the SSBO-approach above, I just provide dummy functions for
    // the BaseStore. So this design didn't age well...
    static void replace(glm::ivec3 group_coord,
                        RaycastEntry* entry)
    {
        entry->world_id = bogus_world_id; // Really didn't age well...
        entry->group_coord = group_coord;
        entry->should_draw = false;
    }

    // Another dummy function.
    static bool update(glm::ivec3, RaycastEntry*, bool)
    {
        return true;
    }

  private:
    // Mark all write staging buffers as unused.
    void mark_write_staging_buffers_unused()
    {
        RaycastStore& store = raycast_store;
        for (StagingBuffer& sb : store.write_staging_buffers) {
            sb.world_id = 0;
        }
    }

    // Handle step 2 of the data flow: try to get a staging buffer
    // (return true iff this is done) and fill it with data from the
    // chunk group with the given group_coordinate.
    //
    // We can clear the dirty flags exactly when this is done
    // successfully.
    bool try_upload_to_staging(glm::ivec3 group_coord)
    {
        RaycastStore& store = raycast_store;
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
        sb.world_id = bogus_world_id;
        sb.group_coord = group_coord;
        sb.make_ready();

        // Load the needed chunk group. No need to upload (report success)
        // if the chunk group does not exist.
        set_current_chunk_group(group_coord);
        if (current_chunk_group == nullptr) return true;

        // Fill in the texels and AABBs. Clear the dirty flag before
        // reading the actual data to handle data races better.
        current_chunk_group->dirty_flags &= ~renderer_raycast_dirty_flag;
        assert(sb.mapped_ssbo != nullptr);
        auto& big_ol_array = sb.mapped_ssbo->voxel_colors;
        StagingBuffer::AABBs a;

        for (int zH = 0; zH < edge_chunks; ++zH) {
        for (int yH = 0; yH < edge_chunks; ++yH) {
        for (int xH = 0; xH < edge_chunks; ++xH) {
            const BinChunk& chunk = current_chunk_group
                                  ->chunk_array[zH][yH][xH];
            glm::ivec3 chunk_residue(
                xH * chunk_size, yH * chunk_size, zH * chunk_size);
            a.aabb_array[zH][yH][xH] = PackedAABB(chunk, chunk_residue);

            for (int zL = 0; zL < chunk_size; ++zL) {
            for (int yL = 0; yL < chunk_size; ++yL) {
            for (int xL = 0; xL < chunk_size; ++xL) {
                big_ol_array[zH][yH][xH][zL][yL][xL] =
                    chunk.voxel_array[zL][yL][xL];
            }
            }
            }
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
        RaycastStore& store = raycast_store;

        auto it = store.queue.begin();

        for ( ; it != store.queue.end(); ++it) {
            auto world_id = it->world_id;
            auto group_coord = it->group_coord;

            if (world_id != bogus_world_id) continue;
            bool success = try_upload_to_staging(group_coord);
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
        RaycastStore& store = raycast_store;
        for (StagingBuffer& sb : store.read_staging_buffers) {
            if (sb.world_id == 0) continue;

            StoreRequest request;
            RaycastEntry* entry =
                store.request_no_update(sb.group_coord, &request);

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
        // glMemoryBarrier(GL_CLIENT_MAPPED_BUFFER_BARRIER_BIT);

        RaycastStore& store = raycast_store;
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
            glDispatchCompute(1, 2, 2);
            PANIC_IF_GL_ERROR;
        }
    }

    // Render, to the current framebuffer, chunks around the camera
    // using the AABB-raycast algorithm. Chunks that are near
    // enough to have been drawn using the mesh algorithm will
    // not be re-drawn.
    void render_world_raycast_step() noexcept
    {
        RaycastStore& store = raycast_store;
        store.stats_begin_frame();

        // Resize the cache to a size suitable for the given render distance.
        auto far_plane = tr.far_plane;
        uint32_t modulus = uint32_t(std::max(4.0, ceil(far_plane / 128.0)));
        store.set_modulus(modulus);
        store.resize_victim_cache(0);

        // Count and limit the number of new chunk groups added to the
        // RaycastStore this frame.
        int remaining_new_chunk_groups_allowed =
            tr.max_frame_new_chunk_groups;

        // Collect all chunk groups needing to be raycast, and sort from
        // nearest to furthest (reverse painters). This fascilitates
        // the early depth test optimization.
        float squared_thresh = tr.far_plane;
        squared_thresh *= squared_thresh;
        std::vector<std::pair<float, glm::ivec3>> group_coord_by_depth;

        unsigned drawn_group_count = 0;
        glm::ivec3 group_coord_low, group_coord_high;
        glm::dvec3 disp(tr.far_plane);
        split_coordinate(tr.eye - disp, &group_coord_low);
        split_coordinate(tr.eye + disp, &group_coord_high);

        for (int32_t zH = group_coord_low.z; zH <= group_coord_high.z; ++zH) {
        for (int32_t yH = group_coord_low.y; yH <= group_coord_high.y; ++yH) {
        for (int32_t xH = group_coord_low.x; xH <= group_coord_high.x; ++xH) {
            auto group_coord = glm::ivec3(xH, yH, zH);
            set_current_chunk_group(group_coord);
            if (current_chunk_group == nullptr) continue;

            float min_squared_dist;
            bool cull = cull_group(group_coord, &min_squared_dist);
            if (cull or min_squared_dist >= squared_thresh) continue;
            group_coord_by_depth.emplace_back(min_squared_dist, group_coord);
            ++drawn_group_count;
        }
        }
        }

        if (!disable_zcull_sort) {
            auto lt_depth = [] (const auto& left, const auto& right)
            {
                return left.first < right.first;
            };
            std::sort(group_coord_by_depth.begin(),
                      group_coord_by_depth.end(),
                      lt_depth);
        }
        std::vector<RaycastEntry*> draw_in_first_pass;
        std::vector<RaycastEntry*> draw_in_second_pass;

        // Implement the earlier pseudocode here.
        store.write_staging_buffers.resize(96);
        store.read_staging_buffers.resize(96);
        mark_write_staging_buffers_unused();
        queue_to_staging_buffers(); // (2)

        for (auto& pair : group_coord_by_depth) {
            // Steps 1 and 2: prepare to add dirty chunk groups to staging.
            auto group_coord = pair.second;
            StoreRequest request;
            request.may_fail = remaining_new_chunk_groups_allowed <= 0;
            RaycastEntry* entry = store.request<Renderer>(
                *this, group_coord, &request);

            if (entry == nullptr) continue;

            set_current_chunk_group(group_coord);
            if (current_chunk_group == nullptr) continue;

            bool dirty = current_chunk_group->dirty_flags
                       & renderer_raycast_dirty_flag;
            if (dirty or !request.was_in_store) {
                bool success = try_upload_to_staging(group_coord); // (2)
                if (!success) {
                    // (1)
                    store.queue.push_back( { bogus_world_id, group_coord } );
                }
            }
            remaining_new_chunk_groups_allowed -= !request.was_in_store;
        }

        swap_read_staging_buffers_into_RaycastStore(); // (6)
        for (auto& pair : group_coord_by_depth) {
            StoreRequest request;
            request.may_fail = true;
            auto group_coord = pair.second;
            RaycastEntry* entry = store.request<Renderer>(
                *this, group_coord, &request);

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
        glm::mat4 residue_vp_matrix = tr.residue_vp_matrix;
        glm::vec3 eye_residue = tr.eye_residue;
        glm::ivec3 eye_group = tr.eye_group;
        auto far_plane = tr.far_plane;

        static GLuint vao = 0;
        static GLuint program_id;
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
        }

        glBindVertexArray(vao);
        glUseProgram(program_id);

        // Set uniforms unchanged per-chunk-group.
        glUniform1i(chunk_debug_id, chunk_debug);
        glUniform1i(fog_enabled_id, tr.use_fog);
        glUniform1i(black_fog_id, tr.use_black_fog);
        glUniform1i(far_plane_squared_id, far_plane * far_plane);
        auto raycast_thr = raycast_threshold;
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

            // Get the instanced vertex attribs going (i.e. bind the
            // AABB residue coords, packed as integers).
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
            // As specified in raycast.vert, need 14 triangle strip drawn.
            auto instances = edge_chunks * edge_chunks * edge_chunks;
            glDrawArraysInstanced(
                GL_TRIANGLE_STRIP, 0, 14, instances);
            PANIC_IF_GL_ERROR;
        }

        PANIC_IF_GL_ERROR;
        glBindVertexArray(0);
    }

    // Render the background. This uses a shader hard-wired to draw a
    // full-screen rectangle.
    void render_background()
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

        glm::vec3 eye_residue = tr.eye_residue;
        glm::mat4 inverse_vp = glm::inverse(tr.residue_vp_matrix);

        glBindVertexArray(vao);
        glUseProgram(program_id);
        glUniformMatrix4fv(inverse_vp_id, 1, 0, &inverse_vp[0][0]);
        glUniform3fv(eye_world_position_id, 1, &eye_residue[0]);
        glUniform1i(black_fog_id, tr.use_black_fog);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glBindVertexArray(0);
        PANIC_IF_GL_ERROR;
    };

    void draw_frame()
    {
        tr = CameraTransforms(*sync_camera_ptr);
        int x = tr.frame_x, y = tr.frame_y;

        glViewport(0, 0, x, y);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
        if (tr.target_fragments > 0) {
            fprintf(stderr, "TODO thread safety.\n");
            bind_global_f32_depth_framebuffer(x, y, tr.target_fragments);
        }
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
        render_world_mesh_step();
        render_world_raycast_step();
        render_background();
        if (tr.target_fragments > 0) {
            finish_global_f32_depth_framebuffer(x, y);
        }

        window_ptr->swap_buffers();
    }

    static void render_loop(Renderer* renderer)
    {
        renderer->window_ptr->gl_make_current();

        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glProvokingVertex(GL_FIRST_VERTEX_CONVENTION);

        // In principle we just have to draw frames in a loop now, but
        // this function is bigger than you expect since I want to
        // compute FPS.

        // Seconds (since glfw initialization?) of previous drawn frame.
        double previous_update = glfwGetTime();

        double previous_fps_update = glfwGetTime();
        int frames = 0;
        double frame_time = 0;

        while (!renderer->thread_exit_flag.load()) {
            // First, actually draw the frame.
            renderer->draw_frame();

            // Calculate (approximate) frame time. Require at least 2
            // ms between frames (workaround to system freeze bug).
            double now;
            double dt;
            do {
                now = glfwGetTime();
                dt = now - previous_update;
            } while (dt < 0.002);
            previous_update = now;

            // Update FPS and frame time. Periodically report the
            // frame time back to the Renderer.
            ++frames;
            frame_time = std::max(frame_time, dt);
            if (now - previous_fps_update >= 0.5) {
                renderer->fps.store(frames / (now - previous_fps_update));
                renderer->frame_time.store(frame_time);

                previous_fps_update = now;
                frames = 0;
                frame_time = 0;
            }
        }
    }
};

Renderer* new_renderer(
    std::shared_ptr<Window> window,
    WorldHandle world,
    std::shared_ptr<SyncCamera> camera)
{
    return new Renderer(std::move(window), std::move(world), std::move(camera));
}

void delete_renderer(Renderer* renderer)
{
    delete renderer;
}

double get_fps(const Renderer& renderer)
{
    return renderer.fps.load();
}

double get_frame_time(const Renderer& renderer)
{
    return renderer.frame_time.load();
}

static GLuint f32_depth_framebuffer = 0;
static GLuint f32_depth_texture = 0;
static GLuint f32_depth_renderbuffer = 0;
static int f32_depth_framebuffer_x, f32_depth_framebuffer_y = 0;

// Used to detect requested screen (framebuffer) size changes.
static int prev_frame_x = 0, prev_frame_y = 0, prev_target_fragments = 0;

// Ratio of screen size (length / width) to framebuffer size.
static int framebuffer_divisor = 1;


void bind_global_f32_depth_framebuffer(
    int frame_x, int frame_y, int target_fragments)
{
    int want_framebuffer_x = f32_depth_framebuffer_x;
    int want_framebuffer_y = f32_depth_framebuffer_y;

    if (target_fragments <= 0) {
        want_framebuffer_x = frame_x;
        want_framebuffer_y = frame_y;
        framebuffer_divisor = 1;
    }

    // Recalculate the actual framebuffer size if needed.
    else if (frame_x != prev_frame_x
         or frame_y != prev_frame_y
         or prev_target_fragments != target_fragments) {

        for (framebuffer_divisor = 1;
             framebuffer_divisor <= 16;
             ++framebuffer_divisor) {

            want_framebuffer_x = frame_x / framebuffer_divisor;
            want_framebuffer_y = frame_y / framebuffer_divisor;
            if (want_framebuffer_x * want_framebuffer_y <= target_fragments) {
                break;
            }
        }
    }

    prev_frame_x = frame_x;
    prev_frame_y = frame_y;
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

void finish_global_f32_depth_framebuffer(int frame_x, int frame_y)
{
    assert(f32_depth_framebuffer != 0);
    glBlitNamedFramebuffer(
        f32_depth_framebuffer,
        0, // Write to window framebuffer
        0, 0,
        f32_depth_framebuffer_x, f32_depth_framebuffer_y,
        0, 0,
        frame_x, frame_y,
        GL_COLOR_BUFFER_BIT,
        GL_NEAREST);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    PANIC_IF_GL_ERROR;
}

} // end namespace
