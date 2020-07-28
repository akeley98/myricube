






























in vec3 color;
in vec3 residue_coord;
in vec2 uv;

out vec4 out_color;

uniform int far_plane_squared;
uniform vec3 eye_relative_group_origin;

void main() {
    const float d = BORDER_WIDTH;
    bool u_border = (uv.x - floor(uv.x + d) < d);
    bool v_border = (uv.y - floor(uv.y + d) < d);
    float border_fade = (u_border || v_border) ? 0.5 : 1.0;

    vec3 disp = residue_coord - eye_relative_group_origin;
    float dist_squared = dot(disp, disp);
    float raw_fog_fade = FOG_SCALAR * (1 - dist_squared/far_plane_squared);
    float fog_fade = clamp(raw_fog_fade, 0, 1);

    out_color = vec4(color * border_fade * fog_fade, 1);
}

