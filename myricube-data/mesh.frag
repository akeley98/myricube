





























in vec3 v_color;
in vec3 v_residue_coord;
out vec4 out_color;
uniform int far_plane_squared;
uniform vec3 eye_relative_group_origin;
uniform bool fog_enabled;

void main() {
    const float d = BORDER_WIDTH;
    float x = v_residue_coord.x;
    int x_border = (x - floor(x + d) < d) ? 1 : 0;
    float y = v_residue_coord.y;
    int y_border = (y - floor(y + d) < d) ? 1 : 0;
    float z = v_residue_coord.z;
    int z_border = (z - floor(z + d) < d) ? 1 : 0;
    float border_fade = (x_border + y_border + z_border >= 2) ? 0.5 : 1.0;

    vec3 disp = v_residue_coord - eye_relative_group_origin;
    float dist_squared = dot(disp, disp);
    float raw_fog_fade =
        fog_enabled ?
        FOG_SCALAR * (1 - sqrt(dist_squared/far_plane_squared)) :
        1.0;
    float fog_fade = clamp(raw_fog_fade, 0, 1);

    out_color = vec4(v_color * border_fade * fog_fade, 1);
}
