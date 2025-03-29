#ifndef MYRICUBE_DITHER_GLSL_
#define MYRICUBE_DITHER_GLSL_

#include "srgb.glsl"

// https://github.com/nvpro-samples/vk_timeline_semaphore/blob/main/shaders/background.frag
// http://www.thetenthplanet.de/archives/5367
// Apply dithering to hide banding artifacts.
// Returns rgb vector encoded as srgb u8 (0-255).
uvec3 srgb8_dither(vec3 linear_color)
{
    // Copied pseudo random number generation code.
    // http://www.jcgt.org/published/0009/03/02/
    // Hash Functions for GPU Rendering, Mark Jarzynski, Marc Olano, NVIDIA
    uvec3 v = uvec3(uvec2(gl_FragCoord.xy), 20010106);
    v = v*1664525u + uvec3(1013904223u);
    v.x += v.y*v.z; v.y += v.z*v.x; v.z += v.x*v.y;
    v ^= v >> uvec3(16u);
    v.x += v.y*v.z; v.y += v.z*v.x; v.z += v.x*v.y;

    vec3 noise = (1.0 / 4294967296.0) * vec3(v);
    uvec3 lowQuant;
    lowQuant.r = srgb8_from_linear_bias(linear_color.r, 0.0);
    lowQuant.g = srgb8_from_linear_bias(linear_color.g, 0.0);
    lowQuant.b = srgb8_from_linear_bias(linear_color.b, 0.0);
    uvec3 highQuant  = lowQuant + uvec3(1);
    vec3  lowLinear  = linear_from_srgb8_vec(lowQuant);
    vec3  highLinear = linear_from_srgb8_vec(highQuant);
    vec3  discr      = mix(lowLinear, highLinear, noise);
    return lowQuant + uvec3(lessThan(discr, linear_color));
    // return mix(lowLinear, highLinear, lessThan(discr, linear_color));
}

#endif