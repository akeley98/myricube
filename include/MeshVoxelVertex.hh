// Utility for converting chunks into a mesh to draw.  I'm using
// instanced rendering, so every drawn voxel will be represented as a
// packed bitfield of its position, color, and which faces are
// visible, instead of as a raw triangles mesh.

#ifndef MYRICUBE_MESHVOXELVERTEX_HH_
#define MYRICUBE_MESHVOXELVERTEX_HH_

#include "myricube.hh"
#include "PackedAABB.hh"
#include "voxels.hh"

namespace myricube {

constexpr uint32_t pos_x_face_bit = (1 << 24);
constexpr uint32_t neg_x_face_bit = (1 << 25);
constexpr uint32_t pos_y_face_bit = (1 << 26);
constexpr uint32_t neg_y_face_bit = (1 << 27);
constexpr uint32_t pos_z_face_bit = (1 << 28);
constexpr uint32_t neg_z_face_bit = (1 << 29);
constexpr uint32_t all_face_bits = pos_x_face_bit
                                 | neg_x_face_bit
                                 | pos_y_face_bit
                                 | neg_y_face_bit
                                 | pos_z_face_bit
                                 | neg_z_face_bit;

// Packed-bitfield "Vertex" for one voxel in a mesh.
struct MeshVoxelVertex
{
    // Packed bitfield of x/y/z residue coordinates, and which of the
    // +/- x/y/z faces are visible.
    uint32_t packed_residue_face_bits = 0xFFFFFFFF;
    // Packed 8-bit red/green/blue.
    uint32_t packed_color = 0xFFFFFFFF;

    MeshVoxelVertex() = default;

    MeshVoxelVertex(uint32_t a, uint32_t b)
    {
        packed_residue_face_bits = a;
        packed_color = b;
    }

    // Given a _visible_ voxel and its residue coordinate, return the
    // VBO vertex carrying this information.  None of the face bits
    // are set (so this voxel starts as invisible).
    MeshVoxelVertex(uint32_t packed_color_arg, uint8_t x, uint8_t y, uint8_t z)
    {
        static_assert(group_size <= 255,
                      "group too big for 8-bit unsigned coordinates.");
        packed_residue_face_bits = uint32_t(x) << x_shift
                                 | uint32_t(y) << y_shift
                                 | uint32_t(z) << z_shift;
        assert(packed_color_arg & visible_bit);
        packed_color = packed_color_arg;
    }
};

// The maximum number of MeshVoxelVerts needed for one chunk group.
// it's just the total number of voxels per chunk. I'm sure there's a
// lower possible bound but for now I'll be conservative (even though
// I'm desperate for GPU memory).
constexpr size_t group_max_verts = group_size * group_size * group_size;

// Bytes on GPU for storing the mesh (list of visible voxels) of one
// chunk group.
struct MappedGroupMesh
{
    MeshVoxelVertex verts[group_max_verts];
};

// Extra data needed to interpret MappedChunkMesh correctly for drawing.
struct ChunkDrawData
{
    // Vertex (visible voxel) counts( i.e. # instances)
    // stored in MappedGroupMesh::verts[first_voxel : first_voxel+voxel_count]
    uint32_t first_voxel = 0, voxel_count = 0;

    PackedAABB aabb;       // AABB, for decide_chunk's benefit.
};

// Function for filling the above structures given a chunk.
// Returns the number of voxels filled in MappedGroupMesh
inline size_t fill_chunk_group_mesh(
    MappedGroupMesh* mesh_ptr,
    ChunkDrawData (&chunk_draw_data)[edge_chunks][edge_chunks][edge_chunks],  // [z][y][x]
    const BinChunkGroup& chunk_group)
{
    uint32_t total_voxels = 0;

    // Look up whether the voxel at the given coordinate
    // (relative to the lower-left of this chunk group) is visible.
    // Act as if voxels outside the chunk group are always invisible.
    auto visible_block = [&chunk_group] (glm::ivec3 coord) -> bool
    {
        if (coord.x < 0 or coord.x >= group_size
         or coord.y < 0 or coord.y >= group_size
         or coord.z < 0 or coord.z >= group_size) return false;
        return chunk_group(coord) & visible_bit;
    };

    auto visit_voxel = [mesh_ptr, &chunk_group, visible_block, &total_voxels]
    (glm::ivec3 coord, ChunkDrawData* draw_data_ptr)
    {
        auto v = chunk_group(coord);
        if (0 == (v & visible_bit)) return;

        uint8_t x = uint8_t(coord.x);
        uint8_t y = uint8_t(coord.y);
        uint8_t z = uint8_t(coord.z);

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
            const auto idx = draw_data_ptr->first_voxel
                             + draw_data_ptr->voxel_count++;
            assert(idx < group_max_verts);
            mesh_ptr->verts[idx] = vert;
            assert(idx == total_voxels);
            total_voxels++;
        }
    };

    for (int zC = 0; zC < edge_chunks; ++zC) {
        for (int yC = 0; yC < edge_chunks; ++yC) {
            for (int xC = 0; xC < edge_chunks; ++xC) {
                glm::ivec3 chunk_index(xC, yC, zC);
                ChunkDrawData* draw_data_ptr = &chunk_draw_data[zC][yC][xC];
                BinChunkView chunk_view{&chunk_group, chunk_index};
                draw_data_ptr->first_voxel = total_voxels;
                draw_data_ptr->voxel_count = 0;
                draw_data_ptr->aabb = PackedAABB(chunk_view);
                for (int zV = 0; zV < chunk_size; ++zV) {
                    for (int yV = 0; yV < chunk_size; ++yV) {
                        for (int xV = 0; xV < chunk_size; ++xV) {
                            auto coord = glm::ivec3(xV, yV, zV)
                                         + chunk_index * chunk_size;
                            visit_voxel(coord, draw_data_ptr);
                        }
                    }
                }
            }
        }
    }
    return total_voxels;
}

} // end namespace

#endif
