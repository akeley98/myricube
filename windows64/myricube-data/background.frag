






























uniform mat4 inverse_vp;
uniform vec3 eye_world_position;

in vec2 normalized_screen_xy;
out vec4 color;

vec3 fog_color_from_world_direction(vec3);

// Basically, reverse engineer this fragment's position and compare it
// to the eye position to get a direction-we're-looking at, which can
// be converted to a fog color.
void main() {
    color = vec4(1.0, 0.0, 1.0, 1.0);
    vec4 v = inverse_vp * vec4(normalized_screen_xy, 1.0, 1.0);
    vec3 frag_world_position = v.xyz / v.w;
    vec3 direction = frag_world_position - eye_world_position;
    vec3 fog_color = fog_color_from_world_direction(direction);
    color = vec4(fog_color, 1.0);
}

