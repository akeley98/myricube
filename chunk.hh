// CPU-side implementations of voxel grid data structures: voxels
// themselves, chunks, chunk groups, and entire "infinite" voxel
// worlds.

#ifndef MYRICUBE_CHUNK_HH_
#define MYRICUBE_CHUNK_HH_

#include "myricube.hh"

#include <assert.h>
#include <atomic>
#include <unordered_map>

namespace myricube {

class Renderer;
class ChunkGroup;

// A texel of the 3D texture used to represent the voxel grid on the GPU.
//
// This will probably become 8-bit palettized color if I really start
// to optimize aggressively.
struct VoxelTexel
{
    uint8_t red = 0, green = 0, blue = 0, alpha = 0;
};

// CPU-side Voxel storage. Just a color for now, but more may come
// later (hence the distinction from the VoxelTexel struct).
struct Voxel
{
    bool visible = false;
    uint8_t red = 0, green = 0, blue = 0;
    operator bool() const { return visible; }

    Voxel() = default;
    Voxel(uint8_t red_, uint8_t green_, uint8_t blue_)
    {
        visible = true;
        red = red_;
        green = green_;
        blue = blue_;
    }

    explicit operator VoxelTexel() const
    {
        VoxelTexel t;
        t.red = red;
        t.green = green;
        t.blue = blue;
        t.alpha = visible ? 255 : 0;
        return t;
    }
};

// Small cubic collection of voxels, with its associated (copy of) the
// 3D texture representing it on the graphics card.
class Chunk
{
    friend class Renderer;
    friend class ChunkGroup;

    // True when the AABB needs to be recalculated.
    bool aabb_dirty = true;

    // True when the AABB needs to be re-sent to the GPU.
    bool aabb_gpu_dirty = true;

    // True when the texture needs to be re-uploaded to the GPU.
    bool texture_dirty = true;

    // True when the mesh needs to be re-uploaded to the GPU.
    bool mesh_dirty = true;

    // As an optimization, consider moving the dirty flags to an array
    // of flags in ChunkGroup, shared among all Chunks in the group.
    // This would improve memory locality.

    // Array of voxels within the chunk, in [z][y][x] order (to match
    // GPU texture).
    Voxel voxel_array[chunk_size][chunk_size][chunk_size];

    // Mirror image of voxel_array, but with texels instead.
    VoxelTexel voxel_texture[chunk_size][chunk_size][chunk_size];

    // {xyz}_visible[n] is the number of visible voxels in this chunk
    // with {xyz} == n.
    int16_t x_visible[chunk_size] = { 0 };
    int16_t y_visible[chunk_size] = { 0 };
    int16_t z_visible[chunk_size] = { 0 };

    // Total number of visible voxels in this chunk.
    int16_t total_visible = 0;
    static_assert(chunk_size * chunk_size * chunk_size <= 32767,
        "int16_t may not be enough for Chunk::total_visible.");

    // residue coordinate of the voxel at voxel_array[0][0][0], i.e.
    // the offset in voxel units that the lower-left corner of this
    // chunk is relative to the lower-left corner of the chunk group
    // that it is in. This needs to be set by ChunkGroup.
    glm::ivec3 chunk_offset;

    // Lower-left corner and upper-right of the AABB in residue
    // coordinates.
    glm::ivec3 aabb_low;
    glm::ivec3 aabb_high;

    // Assuming that the given coordinate is wihin this chunk, return
    // the voxel at the given coordinate. ("Assume" ==
    // the upper bits are masked away without checking).
    Voxel operator() (glm::ivec3 c) const
    {
        auto mask = chunk_size - 1;
        return voxel_array[c.z & mask][c.y & mask][c.x & mask];
    }

    // Set the voxel at the given coordinate (again assuming that this
    // voxel is in the chunk). Return the change in the number of
    // visible voxels in the chunk.
    int set(glm::ivec3 c, Voxel v)
    {
        aabb_dirty = true;
        aabb_gpu_dirty = true;
        texture_dirty = true;
        mesh_dirty = true;

        auto mask = chunk_size-1;
        auto x = mask & c.x;
        auto y = mask & c.y;
        auto z = mask & c.z;

        // Change in the number of visible voxels in the chunk.
        // Use this to update the visible voxel counts.
        int16_t delta = int16_t(v.visible)
                      - int16_t(voxel_array[z][y][x].visible);
        x_visible[x] += delta;
        y_visible[y] += delta;
        z_visible[z] += delta;
        total_visible += delta;
        voxel_array[z][y][x] = v;
        voxel_texture[z][y][x] = VoxelTexel(v);
        return delta;
    }

    // Write out the AABB lower-left corner and upper-right corner,
    // recomputing if needed.
    //
    // This (and the voxel_texture) really belongs in the Renderer, maybe
    // I'll move it some time.
    void get_aabb(glm::ivec3* ptr_aabb_low, glm::ivec3* ptr_aabb_high)
    {
        if (aabb_dirty) {
            if (total_visible == 0) {
                aabb_low = glm::ivec3(0,0,0);
                aabb_high = glm::ivec3(0,0,0);
            }
            else {
                aabb_low = glm::ivec3(chunk_size, chunk_size, chunk_size);
                aabb_high = glm::ivec3(0, 0, 0);

                for (int n = 0; n < chunk_size; ++n) {
                    if (x_visible[n] > 0) {
                        aabb_low.x = std::min<int>(n, aabb_low.x);
                        aabb_high.x = std::max<int>(n+1, aabb_high.x);
                    }
                    if (y_visible[n] > 0) {
                        aabb_low.y = std::min<int>(n, aabb_low.y);
                        aabb_high.y = std::max<int>(n+1, aabb_high.y);
                    }
                    if (z_visible[n] > 0) {
                        aabb_low.z = std::min<int>(n, aabb_low.z);
                        aabb_high.z = std::max<int>(n+1, aabb_high.z);
                    }
                }
            }
            aabb_dirty = false;
            aabb_low += chunk_offset;
            aabb_high += chunk_offset;
        }
        if (ptr_aabb_low) *ptr_aabb_low = aabb_low;
        if (ptr_aabb_high) *ptr_aabb_high = aabb_high;
    }
};

// Cubic collection of chunks. There's not much really going on here
// on the CPU side; this is mostly for rendering purposes.
class ChunkGroup
{
    friend class Renderer;

    // Chunks within this chunk group, in [z][y][x] order.
    Chunk chunk_array[edge_chunks][edge_chunks][edge_chunks];

    // Total number of visible voxels in this chunk group.
    int32_t total_visible = 0;

  public:
    ChunkGroup()
    {
        for (int z = 0; z < edge_chunks; ++z) {
            for (int y = 0; y < edge_chunks; ++y) {
                for (int x = 0; x < edge_chunks; ++x) {
                    chunk_array[z][y][x].chunk_offset =
                        chunk_size * glm::ivec3(x, y, z);
                }
            }
        }
    }

    // Assuming that the given coordinate is wihin this group, return
    // the voxel at the given coordinate.
    Voxel operator() (glm::ivec3 c) const
    {
        auto mask = uint32_t(edge_chunks-1);
        auto x = (uint32_t(c.x) / chunk_size) & mask;
        auto y = (uint32_t(c.y) / chunk_size) & mask;
        auto z = (uint32_t(c.z) / chunk_size) & mask;
        return (chunk_array[z][y][x])(c);
    }

    // Set the voxel at the given coordinate (again assuming that
    // this voxel is in the group).
    void set(glm::ivec3 c, Voxel v)
    {
        auto mask = uint32_t(edge_chunks-1);
        auto x = (uint32_t(c.x) / chunk_size) & mask;
        auto y = (uint32_t(c.y) / chunk_size) & mask;
        auto z = (uint32_t(c.z) / chunk_size) & mask;
        Chunk& chunk = chunk_array[z][y][x];
        total_visible += chunk.set(c, v);
    }
};

// Chunk group labeled with its group coordinate.
//
// Group coordinate = floor(c / group_size), where c is the coordinate
// of any voxel in the group.
using PositionedChunkGroup = std::pair<const glm::ivec3, ChunkGroup>;

// Given a voxel coordinate, convert it to the coordinate of the group
// it is within.
inline glm::ivec3 group_coord(glm::ivec3 voxel_coord)
{
    // Requires arithmetic shift.
    return glm::ivec3(voxel_coord.x >> group_shift,
                      voxel_coord.y >> group_shift,
                      voxel_coord.z >> group_shift);
}

// Convenience / understandability getters for positioned chunk groups.
inline glm::ivec3 group_coord(const PositionedChunkGroup& g)
{
    return g.first;
}

inline ChunkGroup& group(PositionedChunkGroup& g)
{
    return g.second;
}

inline const ChunkGroup& group(const PositionedChunkGroup& g)
{
    return g.second;
}

// Get the voxel with the given coordinate, assuming the coordinate
// is in the group referenced.
inline Voxel get(const PositionedChunkGroup& g, glm::ivec3 c)
{
    assert(group_coord(c) == group_coord(g));
    return group(g)(c);
}

// Set the voxel with the given coordinate, assuming the coordinate
// is in the group referenced.
inline void write(PositionedChunkGroup& g, glm::ivec3 c, Voxel v)
{
    assert(group_coord(c) == group_coord(g));
    group(g).set(c, v);
}

// Entire world of voxels, built out of chunk groups.
class VoxelWorld
{
    friend class Renderer;

    // Stupid hasher for ivec3.
    struct Hash
    {
        size_t operator() (glm::ivec3 v) const
        {
            constexpr size_t mask = (1 << 21) - 1;
            return (size_t(v.z) & mask) |
                   (size_t(v.y) & mask) << 21 |
                   size_t(v.z) << 42;
        }
    };

    // ID of this world, unique within a single run of this program.
    uint64_t unique_id;

    // Chunk groups mapped by group coordinate.
    //
    // I believe that hint pointers to PositionedChunkGroups within
    // are not invalidated as long as the position key is not removed.
    std::unordered_map<glm::ivec3, ChunkGroup, Hash> group_map;

  public:
    VoxelWorld()
    {
        static std::atomic<uint64_t> next_id;
        unique_id = next_id++;
    }

    VoxelWorld(VoxelWorld&&) = delete;

    uint64_t id() const
    {
        return unique_id;
    }

    // Return voxel at the given coordinate. If a pointer to a
    // non-null hint pointer is provided, the chunk pointed to by the
    // hint is checked first. If this hint was incorrect, the correct
    // hint (if any) is written out to *p_hint (if not null).
    //
    // There's a lot of ifs but many will either be inlined out (null checks)
    // or be very predictable in the common case when the hint is correct
    // (as collections of voxels users are interested in tend to be near).
    Voxel operator() (glm::ivec3 c,
                      const PositionedChunkGroup** p_hint = nullptr) const
    {
        glm::ivec3 gc = group_coord(c);
        if (p_hint != nullptr) {
            const PositionedChunkGroup* hint = *p_hint;
            if (hint != nullptr && group_coord(*hint) == gc) {
                return get(*hint, c);
            }
        }
        auto it = group_map.find(gc);
        if (it != group_map.end()) {
            if (p_hint != nullptr) {
                *p_hint = &*it;
            }
            return get(*it, c);
        }
        return Voxel();
    }

    // Same as above, with non-const hint. Needed to cope with C++'s
    // strict const rules (which are justified, hence the removal of
    // the const qualifier here).
    Voxel operator() (glm::ivec3 c,
                      PositionedChunkGroup** p_hint)
    {
        auto cp_hint = const_cast<const PositionedChunkGroup**>(p_hint);
        return (*this)(c, cp_hint);
    }

    // Set the voxel at the given coordinate, with an optional hint as
    // in the above function.
    void set(glm::ivec3 c, Voxel v,
             PositionedChunkGroup** p_hint = nullptr)
    {
        glm::ivec3 gc = group_coord(c);
        if (p_hint != nullptr) {
            PositionedChunkGroup* hint = *p_hint;
            if (hint != nullptr && group_coord(*hint) == gc) {
                write(*hint, c, v);
                return;
            }
        }
        auto it = group_map.try_emplace(gc).first;
        write(*it, c, v);
        if (p_hint != nullptr) *p_hint = &*it;
    }

    // Return voxel at the given coordinate, checking the pointed-to
    // chunk group first. (Maybe make this more robust?)
    Voxel operator() (glm::ivec3 c, const PositionedChunkGroup& hint) const
    {
        const PositionedChunkGroup* p_hint = &hint;
        return (*this)(c, &p_hint);
    }
};

} // end namespace
#endif /* !MYRICUBE_CHUNK_HH_ */
