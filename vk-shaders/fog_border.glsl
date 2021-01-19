#ifndef MYRICUBE_FOG_BORDER_GLSL_
#define MYRICUBE_FOG_BORDER_GLSL_

#include "PushConstant.glsl"

#define BORDER_WIDTH_LOW  0.075
#define BORDER_WIDTH_HIGH 0.45
#define BORDER_FADE  0.75
#define BORDER_DIST_LOW  100
#define BORDER_DIST_HIGH 350
#define DITHER_MASK 3
#define DITHER_SCALAR (0.0035 / DITHER_MASK)
// #define FOG_SCALAR 1.125

float logistic(float t, float k)
{
    return 1.0 / (1 + exp(-k * t));
}

// Convert world-space direction vector into a fog color. (i.e. this
// is the sky color)
vec3 fog_color_from_world_direction(vec3 world_direction)
{
    if ((pc.pc.flags & MYRICUBE_BLACK_FOG_BIT) != 0) {
        return vec3(0, 0, 0);
    }

    world_direction = normalize(world_direction);
    float t_horizon = logistic((world_direction.y + .25), -5.0);
    float t_sunrise = logistic((world_direction.y - 1.5*world_direction.x), 0.8);

    const vec3 earth_color = vec3(0.03, 0.023, 0.074);
    const vec3 sky_color = vec3(0.242, 0.651, 0.939);
    // const vec3 sunrise_color = vec3(1.0, 0.418, 0.254);
    const vec3 sunrise_color = vec3(0.708, 0.445, 0.638);

    const ivec2 pixel = ivec2(gl_FragCoord.xy);
    const vec3 dither = vec3(
        float(pixel.x & DITHER_MASK) * DITHER_SCALAR,
        float((pixel.x ^ pixel.y) & DITHER_MASK) * DITHER_SCALAR,
        float(pixel.y & DITHER_MASK) * DITHER_SCALAR);

    return mix(
        mix(sunrise_color, sky_color, t_sunrise),
        earth_color, t_horizon) + dither;
}

// Utility function for adding fog and border effects (and setting
// alpha=1).  This runs in srgb color space as my original magic
// numbers and equations (that I don't undertand anymore) were tweaked
// in srgb, not linear color. (TODO: I ignore this fact for now)
vec4 fog_border_color(
    vec3 base_color, // The stored color of the voxel.
    float dist_squared, // Squared distance from eye to this fragment.
    vec2 uv, // "Texture coordinate"
    vec3 fog_color)
{
    // Border fade diminishes with distance. First, calculate how
    // strong the border fade is (which might not actually matter if
    // this fragment is not on the border).
    float dist = sqrt(dist_squared);
    const float slope = (1-BORDER_FADE) / (BORDER_DIST_HIGH - BORDER_DIST_LOW);
    float base_border_fade = clamp(
        BORDER_FADE + (dist - BORDER_DIST_LOW) * slope,
        BORDER_FADE, 1.0);

    // Calculate how close this fragment is to the edge of the voxel.
    // Actual border fade is based on this.
    float u_center_dist = abs(uv.x - 0.5 - floor(uv.x));
    float v_center_dist = abs(uv.y - 0.5 - floor(uv.y));
    float l4 = u_center_dist * u_center_dist * u_center_dist * u_center_dist
             + v_center_dist * v_center_dist * v_center_dist * v_center_dist;
    const float magic_low = (0.5-BORDER_WIDTH_LOW) * (0.5-BORDER_WIDTH_LOW)
                          * (0.5-BORDER_WIDTH_LOW) * (0.5-BORDER_WIDTH_LOW);
    const float magic_high = (0.5-BORDER_WIDTH_HIGH) * (0.5-BORDER_WIDTH_HIGH)
                           * (0.5-BORDER_WIDTH_HIGH) * (0.5-BORDER_WIDTH_HIGH);
    float border_fade = clamp(
        1 + (l4 - magic_high) * (base_border_fade - 1)
        * (1 / (magic_low - magic_high)),
        base_border_fade, 1.0);

    float squared_dist_ratio = dist_squared/pc.pc.far_plane_squared;
    float dist_ratio = sqrt(dist_squared/pc.pc.far_plane_squared);
    float raw_fog_fade =
        ((pc.pc.flags & MYRICUBE_FOG_BIT) != 0) ?
        1.01357 * logistic(squared_dist_ratio - 0.5, -10.0) - 0.007 :
        // 1.1 * logistic(squared_dist_ratio - 0.5, -6.0) - 0.05 :
        // 1.1 * logistic(dist_ratio - 0.5, -6.0) - 0.05 :
        // 2.164 * logistic(squared_dist_ratio - 0.5, -2.0) - 0.582 :
        1.0;
    float fog_fade = clamp(raw_fog_fade, 0.0, 1.0);

    // Apply fog and border effects.
    return
    vec4(fog_fade * border_fade * base_color + (1-fog_fade) * fog_color, 1.0);
}

#endif
