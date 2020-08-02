






























#define BORDER_WIDTH 0.25
#define BORDER_FADE  0.6
#define BORDER_DIST_LOW  130
#define BORDER_DIST_HIGH 340
uniform bool fog_enabled;
uniform int far_plane_squared;

// Utility function for adding fog and border effects (and setting alpha=1).
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
        BORDER_FADE, 1);

    // Calculate how close this fragment is to the edge of the voxel.
    // Actual border fade is based on this.
    float u_border_dist = abs(uv.x - floor(uv.x + 0.5));
    float v_border_dist = abs(uv.y - floor(uv.y + 0.5));
    float border_param = clamp(
        1 - min(u_border_dist, v_border_dist) * (1 / BORDER_WIDTH),
        0, 1);
    float border_fade = 1 - border_param * (1-base_border_fade);

    // Fog fade depends on linear distance from camera.
    float raw_fog_fade =
        fog_enabled ?
        FOG_SCALAR * (1 - sqrt(dist_squared/far_plane_squared)) :
        1.0;
    float fog_fade = clamp(raw_fog_fade, 0, 1);

    // Apply fog and border effects.
    return
    vec4(fog_fade * border_fade * base_color + (1-fog_fade) * fog_color, 1);
}
