// Common PushConstant structure shared between host vulkan code and
// shaders. Need some preprocessor hackery to share this.

#ifndef MYRICUBE_PUSHCONSTANT_GLSL_
#define MYRICUBE_PUSHCONSTANT_GLSL_

#ifdef __cplusplus
// Host code defines
#include <stdint.h>
#define MYRICUBE_GLM glm::
#else
// Shader code defines
#define int32_t int
#define MYRICUBE_GLM
#endif

struct PushConstant
{
    MYRICUBE_GLM mat4 mvp;
    MYRICUBE_GLM vec4 eye_relative_group_origin; // w unused.
    int32_t flags;
    int32_t far_plane_squared;
    int32_t raycast_thresh_squared;
};

#define MYRICUBE_FOG_BIT 1
#define MYRICUBE_BLACK_FOG_BIT 2
#define MYRICUBE_CHUNK_DEBUG_BIT 4

#ifdef __cplusplus
static_assert(sizeof(PushConstant) <= 128,
    "Vulkan spec only guarantees 128 bytes of push constant");
#else
layout (push_constant) uniform PushConstantBlock {
    PushConstant pc;
} pc;
#endif

#endif

