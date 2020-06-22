// Common include file for myricube.
//
// Myricube: Experimental voxel renderer using hybrid raycast/mesh
// OpenGL rendering.
//
// David Zhao Akeley 2020 (not a great year)
//
// This was inspired by my discontent with Minecraft's low render
// distance, even on high settings.
//
// Like most voxel renderers, this renderer subdivides an "infinite"
// grid of voxels into cubic chunks. Nearby chunks are drawn using
// conventional methods (by converting the voxels into a mesh of
// triangles). Farther chunks are instead drawn using raycasting: a
// minimal AABB is calculated that contains all the solid voxels of
// the chunk, and the AABB itself is drawn using a raycasting fragment
// shader that checks for collisions with the voxels contained in the
// AABB. The voxel data itself (when raycasting is used) is stored
// using a 3D texture -- this is considerably more memory efficient
// than a triangle mesh.
//
// Cubes of chunks are organized into larger chunk groups. For OpenGL
// efficiency, chunks within the same chunk group share GPU resources
// (3D textures, buffer storage, etc.). Points in space are frequently
// expressed as "group" and "residue" coordinates -- these are the
// coordinates (in chunk-group-size units) of the lower-left of the
// chunk the point is in, and the remaining offset within the
// chunk. (See split_coordinate for the precise floor-based
// definition). For example, if chunk groups were 100 x 100 x 100
// voxels, then point (1, 502.5, -1) has chunk coordinate (0, 5, -1)
// and residue (1, 2.5, 99).
//
// To be precise, note also that a voxel at coordinate (x,y,z) occupies
// the cube from (x,y,z) to (x+1,y+1,z+1) in space.
//
// The primary benefits of raycasting are the aforementioned memory
// efficiency and the reduced number of triangle vertices processed
// (12 triangles are shared for all voxels per chunk, instead of up to
// 12 triangles per solid voxel). There are several drawbacks that make
// raycasting suitable only for distant chunks:
//
// * Raycasting cost is roughly linear in the size of the rendered
//   triangles (because the raycast fragment shader is the most
//   expensive part). This makes raycasting very expensive for large
//   nearby triangles.
//
// * Raycasting is more prone to rounding errors. This creates
//   unsightly gaps between voxels that are more obvious when nearby.
//
// * Raycasting causes "incorrect" writes to the depth buffer (because
//   the depth by default is that of the AABB, and not the actual
//   voxel), so raycast voxels may incorrectly occlude other geometry
//   (imagine a Minecraft mob standing on a voxel in the center of a
//   chunk). There are potential solutions (e.g. manual gl_FragDepth)
//   but the cost in lost optimizations is considerable.
//
// * Circular reasoning, but in general I've designed the raycaster to
//   be cheap-but-ugly.
//
// This motivates the hybrid renderer I've come up with.

#ifndef MYRICUBE_MYRICUBE_HH_
#define MYRICUBE_MYRICUBE_HH_

#define GLM_FORCE_RADIANS
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>

namespace myricube {

// Length (in voxels) of the edge of a chunk group/chunk.
constexpr int group_size = 64;
constexpr int chunk_size = 16;

// Bit shift counts of above.
constexpr int group_shift = 6;
constexpr int chunk_shift = 4;

static_assert((group_size & (group_size-1)) == 0,
    "group_size should be a power of 2.");
static_assert((chunk_size & (chunk_size-1)) == 0,
    "chunk_size should be a power of 2.");
static_assert((group_size % chunk_size) == 0,
    "There must be an integer number of chunks per chunk group.");
static_assert((1 << chunk_shift) == chunk_size, "chunk_shift is wrong");
static_assert((1 << group_shift) == group_size, "group_shift is wrong");

// Edge-length of a chunk group cube in chunks.
constexpr int edge_chunks = group_size / chunk_size;

// Split a 3-vector in voxel units into its group and residue coordinates.
inline void split_coordinate(glm::dvec3 v,
                             glm::ivec3* p_group=nullptr,
                             glm::vec3* p_residue=nullptr)
{
    int32_t group_x = int32_t(glm::floor(v.x / group_size));
    int32_t group_y = int32_t(glm::floor(v.y / group_size));
    int32_t group_z = int32_t(glm::floor(v.z / group_size));

    glm::ivec3 group_coord(group_x, group_y, group_z);
    if (p_group != nullptr) *p_group = group_coord;
    if (p_residue != nullptr) {
        *p_residue =
            glm::vec3(v - glm::dvec3(group_coord) * double(group_size));
    }
}

// Name a file that is in the data directory, and return its absolute path.
std::string expand_filename(const std::string& in);

inline bool is_real(float v)
{
    return v - v == 0;
}

inline bool is_real(double v)
{
    return v - v == 0;
}

inline bool is_real(glm::vec3 v)
{
    return v - v == glm::vec3(0, 0, 0);
}

inline bool is_real(glm::dvec3 v)
{
    return v - v == glm::dvec3(0, 0, 0);
}

inline void panic(const std::string& reason)
{
    fprintf(stderr, "%s\n", reason.c_str());
    exit(1);
}

inline void panic(const char* a, const char* b)
{
    std::string reason = a + std::string("\n") + b;
    panic(reason);
}

#define PANIC_IF_GL_ERROR do { \
    if (GLenum PANIC_error = glGetError()) { \
        char PANIC_msg[160]; \
        snprintf(PANIC_msg, sizeof PANIC_msg, "line %i: code %u", __LINE__, (unsigned)PANIC_error); \
        panic("OpenGL error", PANIC_msg); \
    } \
} while (0)

} // end namespace

#endif /* !MYRICUBE_MYRICUBE_HH_ */
