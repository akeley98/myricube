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

// The maximum number of MeshVoxelVerts needed for one chunk.  For now
// it's just the total number of voxels per chunk. I'm sure there's a
// lower possible bound but for now I'll be conservative (even though
// I'm desperate for GPU memory).
constexpr size_t chunk_max_verts = chunk_size * chunk_size * chunk_size;

// Bytes on GPU for storing the mesh (list of visible voxels) of one
// chunk.
struct MappedChunkMesh
{
    MeshVoxelVertex verts[chunk_max_verts];
};

// An entire chunk group of MappedChunkMesh, [z][y][x] order as typical.
struct MappedGroupMesh
{
    MappedChunkMesh chunks[edge_chunks][edge_chunks][edge_chunks];
};

// Extra data needed to interpret MappedChunkMesh correctly for drawing.
struct ChunkDrawData
{
    size_t vert_count = 0; // Vertex (visible voxel) count i.e. # instances
    PackedAABB aabb;       // AABB, for decide_chunk's benefit.
};

// Function for filling the above structures given a chunk. Requires
// in addition the residue coordinate of the chunk's lower-left corner
// (needed for positioning the AABB and voxel positions in the
// containing chunk group's coordinate system).
inline void fill_chunk_mesh(
    MappedChunkMesh* mesh_ptr,
    ChunkDrawData* draw_data_ptr,
    const BinChunk& chunk,
    glm::ivec3 chunk_residue)
{
    draw_data_ptr->aabb = PackedAABB(chunk, chunk_residue);

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

    auto visit_voxel = [
        mesh_ptr, draw_data_ptr, &chunk, visible_block, chunk_residue]
    (glm::ivec3 coord)
    {
        auto v = chunk(coord);
        if (0 == (v & visible_bit)) return;

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
            assert(draw_data_ptr->vert_count < chunk_max_verts);
            mesh_ptr->verts[draw_data_ptr->vert_count++] = vert;
        }
    };

    draw_data_ptr->vert_count = 0;

    for (int z = 0; z < chunk_size; ++z) {
        for (int y = 0; y < chunk_size; ++y) {
            for (int x = 0; x < chunk_size; ++x) {
                visit_voxel(glm::ivec3(x, y, z));
            }
        }
    }
}

} // end namespace

#endif
