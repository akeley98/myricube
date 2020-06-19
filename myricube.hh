// Common include file for myricube.

#ifndef MYRICUBE_MYRICUBE_HH_
#define MYRICUBE_MYRICUBE_HH_

#define GLM_FORCE_RADIANS
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include <stdint.h>
#include <stdlib.h>

#include <string>

namespace myricube {

constexpr float near_plane = 0.1f, far_plane = 2048.0f;

// Length (in voxels) of the edge of a chunk group/chunk.
constexpr int group_size = 64;
constexpr int chunk_size = 16;

static_assert((group_size & (group_size-1)) == 0,
    "group_size should be a power of 2.");
static_assert((chunk_size & (chunk_size-1)) == 0,
    "chunk_size should be a power of 2.");
static_assert((group_size % chunk_size) == 0,
    "There must be an integer number of chunks per chunk group.");

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

} // end namespace

#endif /* !MYRICUBE_MYRICUBE_HH_ */
