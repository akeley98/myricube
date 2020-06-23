// Common include file for myricube.
//
// See myricube.cc for overview.

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


// Bit shift counts of the below.
constexpr int group_shift = 6;
constexpr int chunk_shift = 4;

// Length (in voxels) of the edge of a chunk group/chunk.
constexpr int group_size = (1 << group_shift);
constexpr int chunk_size = (1 << chunk_shift);

// Edge-length of a chunk group cube in chunks.
constexpr int edge_chunks = group_size / chunk_size;

static_assert((group_size & (group_size-1)) == 0,
    "group_size should be a power of 2.");
static_assert((chunk_size & (chunk_size-1)) == 0,
    "chunk_size should be a power of 2.");
static_assert((group_size % chunk_size) == 0,
    "There must be an integer number of chunks per chunk group.");
static_assert((1 << chunk_shift) == chunk_size, "chunk_shift is wrong");
static_assert((1 << group_shift) == group_size, "group_shift is wrong");


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
