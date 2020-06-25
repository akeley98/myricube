


















layout(location=UNIT_BOX_VERTEX_IDX) in vec4 unit_box_vertex;
layout(location=PACKED_AABB_LOW_IDX) in int packed_aabb_low;
layout(location=PACKED_AABB_HIGH_IDX) in int packed_aabb_high;
out vec3 residue_coord;
out float border_fade;
flat out ivec3 aabb_low;
flat out ivec3 aabb_high;
uniform mat4 mvp_matrix;
uniform vec3 eye_relative_group_origin;
uniform int far_plane_squared;
uniform int raycast_thresh_squared;
void main() {
    int low_x = packed_aabb_low & 255;
    int low_y = (packed_aabb_low >> 8) & 255;
    int low_z = (packed_aabb_low >> 16) & 255;
    int high_x = packed_aabb_high & 255;
    int high_y = (packed_aabb_high >> 8) & 255;
    int high_z = (packed_aabb_high >> 16) & 255;
    vec3 f_aabb_low = vec3(low_x, low_y, low_z);
    vec3 sz = vec3(high_x, high_y, high_z) - f_aabb_low;
    vec4 model_space_pos = vec4(unit_box_vertex.xyz * sz + f_aabb_low, 1);
    vec3 disp = model_space_pos.xyz - eye_relative_group_origin;
    float distance = sqrt(dot(disp, disp));
    border_fade = clamp(distance * 0.003 + 0.11, 0.5, 1.0);
    aabb_low = ivec3(low_x, low_y, low_z);
    aabb_high = ivec3(high_x, high_y, high_z);
    // Re-implement decide_chunk(...) == draw_raycast on GPU.
    vec3 aabb_center = vec3(aabb_low + aabb_high) * 0.5;
    vec3 floor_eye = floor(eye_relative_group_origin);
    disp = aabb_center - floor_eye;
    residue_coord = model_space_pos.xyz;
    float squared_dist = dot(disp, disp);
    bool draw_raycast = squared_dist
        == clamp(squared_dist, raycast_thresh_squared, far_plane_squared);
    // Draw as degenerate triangle if this chunk is not meant for raycasting.
    gl_Position = draw_raycast ? mvp_matrix * model_space_pos
                               : vec4(0,0,0,1);
}
