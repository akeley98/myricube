// Common include file for myricube.
//
// See myricube.cc for overview.

#ifndef MYRICUBE_MYRICUBE_HH_
#define MYRICUBE_MYRICUBE_HH_

#define GLM_FORCE_RADIANS
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include <stdexcept>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>

namespace myricube {

// bit assignments for packed 32-bit colors.
constexpr uint32_t red_shift = 0;
constexpr uint32_t green_shift = 8;
constexpr uint32_t blue_shift = 16;
constexpr uint32_t visible_bit = uint32_t(1) << 31;

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

// Given a voxel coordinate, convert it to the coordinate of the group
// it is within.
inline glm::ivec3 to_group_coord(glm::ivec3 voxel_coord)
{
    // Requires arithmetic shift.
    return glm::ivec3(voxel_coord.x >> group_shift,
                      voxel_coord.y >> group_shift,
                      voxel_coord.z >> group_shift);
}

// Windows uses UTF-16, sigh. (wstring is 16-bit on windows). TODO test.
#if defined(__WIN32__) || defined(__WIN32) || defined(WIN32)
#define MYRICUBE_WINDOWS 1
using filename_string = std::wstring;
#endif

// Meanwhile the free world uses UTF-8 like civilized individuals.
#if defined(__linux__) || defined(__unix__)
using filename_string = std::string;
#endif

using filename_char = filename_string::value_type;

// Dumb function for appending an ascii C string to a filename.
// Needed to deal with platform differences.
inline filename_string filename_concat_c_str(
    filename_string filename, const char* c_str)
{
    auto len = strlen(c_str);
    filename.reserve(filename.size() + len + 1);

    for (size_t i = 0; i < len; ++i) {
        char c = c_str[i];
        // Must be ASCII
        if ((int(c) & 0x7F) != int(c)) {
            throw std::runtime_error(
                "filename_concat_c_str: c_str must be ASCII.\n");
        }
        filename.push_back(static_cast<filename_char>(c));
    }
    return filename;
}

// Name a file that is in the data directory, and return its absolute
// path.  CHANGED BEHAVIOR: no longer detect absolute paths (start
// with /) and leave them unchanged. Arguments now must not be
// absolute paths; they must name a file to be found in the data dir.
filename_string expand_filename(const std::string& in);

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
