





























in vec3 v_color;
in vec3 v_residue_coord;
in vec2 v_uv;

out vec4 out_color;

uniform vec3 eye_relative_group_origin;

vec3 fog_color_from_world_direction(vec3 world_direction);

vec4 fog_border_color(
    vec3 base_color,
    float dist_squared,
    vec2 uv,
    vec3 fog_color);

void main() {
    vec3 disp = v_residue_coord - eye_relative_group_origin;
    float dist_squared = dot(disp, disp);
    vec3 fog_color = fog_color_from_world_direction(disp);
    out_color = fog_border_color(
        v_color, dist_squared, v_uv, fog_color);
}
