#ifndef MYRICUBE_MESHCOMMON_GLSL_
#define MYRICUBE_MESHCOMMON_GLSL_

// Bit assignments for unpacking packed voxels.
#define POS_X_FACE_BIT (1 << 24)
#define NEG_X_FACE_BIT (1 << 25)
#define POS_Y_FACE_BIT (1 << 26)
#define NEG_Y_FACE_BIT (1 << 27)
#define POS_Z_FACE_BIT (1 << 28)
#define NEG_Z_FACE_BIT (1 << 29)

#define X_SHIFT 0
#define Y_SHIFT 8
#define Z_SHIFT 16

#define RED_SHIFT 0
#define GREEN_SHIFT 8
#define BLUE_SHIFT 16

// Other (duplicated) constants.
#define CHUNK_SIZE 16
#define CHUNK_MAX_VOXELS (CHUNK_SIZE * CHUNK_SIZE * CHUNK_SIZE)
#define EDGE_CHUNKS 4

// Needed implementation limits
#define MAX_MESH_VOXELS 32 /* Need NOT be power of 2 */
#define MAX_MESH_VERTICES (MAX_MESH_VOXELS * 8)
#define MAX_MESH_TRIANGLES (MAX_MESH_VOXELS * 12)
#define MAX_TASK_OUTPUTS (CHUNK_MAX_VOXELS / MAX_MESH_VOXELS)

#ifndef __cplusplus

// Instructions for one mesh shader group invocation.
// Which chunk within the chunk group to draw, and offset and count
// within that chunk's verts array to draw.
struct TaskInfo
{
    uint  voxel_offset, voxel_count;
};

// Duplicating definitions in MeshVoxelVertex.hh :/
struct MeshVoxelVertex
{
    uint packed_residue_face_bits;
    uint packed_color;
};

struct ChunkMesh
{
    MeshVoxelVertex verts[CHUNK_MAX_VOXELS];
};

struct ChunkDrawData
{
    uint voxel_count;
    uint packed_low;
    uint packed_high;
};

struct GroupMesh
{
    ChunkMesh        chunks[EDGE_CHUNKS][EDGE_CHUNKS][EDGE_CHUNKS];
    ChunkDrawData draw_data[EDGE_CHUNKS][EDGE_CHUNKS][EDGE_CHUNKS];
};

// Unit cube model, used in shaders only.

// Position coordinate of unit box:
vec3 unit_box_verts[36] = vec3[] (
    vec3(0, 1, 1),    vec3(0, 1, 0),    vec3(0, 0, 1),
    vec3(0, 1, 0),    vec3(0, 0, 0),    vec3(0, 0, 1), // -x face

    vec3(1, 0, 0),    vec3(1, 1, 0),    vec3(1, 0, 1),
    vec3(1, 1, 0),    vec3(1, 1, 1),    vec3(1, 0, 1), // +x face

    vec3(0, 0, 0),    vec3(1, 0, 1),    vec3(0, 0, 1),
    vec3(0, 0, 0),    vec3(1, 0, 0),    vec3(1, 0, 1), // -y face

    vec3(1, 1, 1),    vec3(1, 1, 0),    vec3(0, 1, 0),
    vec3(0, 1, 1),    vec3(1, 1, 1),    vec3(0, 1, 0), // +y face

    vec3(0, 0, 0),    vec3(0, 1, 0),    vec3(1, 0, 0),
    vec3(1, 1, 0),    vec3(1, 0, 0),    vec3(0, 1, 0), // -z face

    vec3(0, 1, 1),    vec3(1, 0, 1),    vec3(1, 1, 1),
    vec3(1, 0, 1),    vec3(0, 1, 1),    vec3(0, 0, 1));// +z face

// Bit corresponding to the face being drawn of the unit box; for
// example, if we're drawing the +x face, this is POS_X_FACE_BIT.
int unit_box_face_bits[36] = int[] (
    NEG_X_FACE_BIT,    NEG_X_FACE_BIT,    NEG_X_FACE_BIT,
    NEG_X_FACE_BIT,    NEG_X_FACE_BIT,    NEG_X_FACE_BIT,

    POS_X_FACE_BIT,    POS_X_FACE_BIT,    POS_X_FACE_BIT,
    POS_X_FACE_BIT,    POS_X_FACE_BIT,    POS_X_FACE_BIT,

    NEG_Y_FACE_BIT,    NEG_Y_FACE_BIT,    NEG_Y_FACE_BIT,
    NEG_Y_FACE_BIT,    NEG_Y_FACE_BIT,    NEG_Y_FACE_BIT,

    POS_Y_FACE_BIT,    POS_Y_FACE_BIT,    POS_Y_FACE_BIT,
    POS_Y_FACE_BIT,    POS_Y_FACE_BIT,    POS_Y_FACE_BIT,

    NEG_Z_FACE_BIT,    NEG_Z_FACE_BIT,    NEG_Z_FACE_BIT,
    NEG_Z_FACE_BIT,    NEG_Z_FACE_BIT,    NEG_Z_FACE_BIT,

    POS_Z_FACE_BIT,    POS_Z_FACE_BIT,    POS_Z_FACE_BIT,
    POS_Z_FACE_BIT,    POS_Z_FACE_BIT,    POS_Z_FACE_BIT);

// Texture coordinate (for now, just used for the color border).
vec2 uv_array[36] = vec2[] (
    vec2(1, 1),    vec2(1, 0),    vec2(0, 1),
    vec2(1, 0),    vec2(0, 0),    vec2(0, 1), // -x face

    vec2(0, 0),    vec2(1, 0),    vec2(0, 1),
    vec2(1, 0),    vec2(1, 1),    vec2(0, 1), // +x face

    vec2(0, 0),    vec2(1, 1),    vec2(0, 1),
    vec2(0, 0),    vec2(1, 0),    vec2(1, 1), // -y face

    vec2(1, 1),    vec2(1, 0),    vec2(0, 0),
    vec2(0, 1),    vec2(1, 1),    vec2(0, 0), // +y face

    vec2(0, 0),    vec2(0, 1),    vec2(1, 0),
    vec2(1, 1),    vec2(1, 0),    vec2(0, 1), // -z face

    vec2(0, 1),    vec2(1, 0),    vec2(1, 1),
    vec2(1, 0),    vec2(0, 1),    vec2(0, 0));// +z face
#endif

#endif
