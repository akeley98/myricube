






























uniform bool fog_enabled;
uniform int far_plane_squared;

// Utility function for adding fog and border effects (and setting alpha=1).
vec4 fog_border_color(
    vec3 base_color, // The stored color of the voxel.
    float dist_squared, // Squared distance from eye to this fragment.
    vec2 uv, // "Texture coordinate"
    float base_border_fade, // Relative brightness of a voxel border.
    vec3 fog_color)
{
    bool u_border = uv.x - floor(uv.x + BORDER_WIDTH) < BORDER_WIDTH;
    bool v_border = uv.y - floor(uv.y + BORDER_WIDTH) < BORDER_WIDTH;
    float border_fade = u_border || v_border ? base_border_fade : 1.0;
    
    float raw_fog_fade =
        fog_enabled ?
        FOG_SCALAR * (1 - sqrt(dist_squared/far_plane_squared)) :
        1.0;
    float fog_fade = clamp(raw_fog_fade, 0, 1);
    
    return
    vec4(fog_fade * border_fade * base_color + (1-fog_fade) * fog_color, 1);
}

