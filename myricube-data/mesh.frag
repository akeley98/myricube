





























in vec3 color;
in vec3 residue_coord;
out vec4 out_color;
uniform int far_plane_squared;
uniform vec3 eye_relative_group_origin;

void main() {
    const float d = BORDER_WIDTH;
    float x = residue_coord.x;
    int x_border = (x - floor(x + d) < d) ? 1 : 0;
    float y = residue_coord.y;
    int y_border = (y - floor(y + d) < d) ? 1 : 0;
    float z = residue_coord.z;
    int z_border = (z - floor(z + d) < d) ? 1 : 0;
    float border_fade = (x_border + y_border + z_border >= 2) ? 0.5 : 1.0;

    vec3 disp = residue_coord - eye_relative_group_origin;
    float dist_squared = dot(disp, disp);
    float raw_fog_fade = FOG_SCALAR * (1 - sqrt(dist_squared/far_plane_squared));
    float fog_fade = clamp(raw_fog_fade, 0, 1);

    out_color = vec4(color * border_fade * fog_fade, 1);
}
