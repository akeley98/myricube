





























// Vertex shader that expands packed AABB (defines lower and upper
// corner of AABB -- passed as instanced vertex attribute) into actual
// AABBs by stretching and positioning a built-in unit box.  Draw as a
// 14-vertex triangle strip, with FIRST vertex as provoking vertex and
// counter-clockwise front (unculled) faces.

// https://stackoverflow.com/questions/28375338/cube-using-single-gl-triangle-strip

vec4 unit_box_verts[14] = vec4[] (
    vec4(1, 1, 1, 1), // +z face
    vec4(0, 1, 1, 1), // +z face
    vec4(1, 0, 1, 1), // -y face part 1
    vec4(0, 0, 1, 1), // -x face
    vec4(0, 0, 0, 1), // -x face
    vec4(0, 1, 1, 1), // +y face
    vec4(0, 1, 0, 1), // +y face
    vec4(1, 1, 1, 1), // +x face
    vec4(1, 1, 0, 1), // +x face
    vec4(1, 0, 1, 1), // -y face part 2
    vec4(1, 0, 0, 1), // -z face
    vec4(0, 0, 0, 1), // -z face
    vec4(1, 1, 0, 1),
    vec4(0, 1, 0, 1));

vec4 unit_box_normals[14] = vec4[] (
    vec4(0, 0, 1, 0), // +z face
    vec4(0, 0, 1, 0), // +z face
    vec4(0,-1, 0, 0), // -y face part 1
    vec4(-1,0, 0, 0), // -x face
    vec4(-1,0, 0, 0), // -x face
    vec4(0, 1, 0, 0), // +y face
    vec4(0, 1, 0, 0), // +y face
    vec4(1, 0, 0, 0), // +x face
    vec4(1, 0, 0, 0), // +x face
    vec4(0,-1, 0, 0), // -y face part 2
    vec4(0, 0,-1, 0), // -z face
    vec4(0, 0,-1, 0), // -z face
    vec4(0, 0, 0, 0),
    vec4(0, 0, 0, 0));

layout(location=PACKED_AABB_LOW_IDX) in int packed_aabb_low;
layout(location=PACKED_AABB_HIGH_IDX) in int packed_aabb_high;
out vec3 aabb_residue_coord;
flat out ivec3 aabb_low;
flat out ivec3 aabb_high;
flat out vec3 aabb_half_normal;
uniform mat4 mvp_matrix;
uniform vec3 eye_relative_group_origin;
uniform int far_plane_squared;
uniform int raycast_thresh_squared;
void main() {
    int low_x = (packed_aabb_low >> X_SHIFT) & 255;
    int low_y = (packed_aabb_low >> Y_SHIFT) & 255;
    int low_z = (packed_aabb_low >> Z_SHIFT) & 255;
    int high_x = (packed_aabb_high >> X_SHIFT) & 255;
    int high_y = (packed_aabb_high >> Y_SHIFT) & 255;
    int high_z = (packed_aabb_high >> Z_SHIFT) & 255;
    vec3 f_aabb_low = vec3(low_x, low_y, low_z);
    vec3 sz = vec3(high_x, high_y, high_z) - f_aabb_low;
    vec3 unit_box_vertex = unit_box_verts[gl_VertexID].xyz;
    vec4 model_space_pos = vec4(unit_box_vertex * sz + f_aabb_low, 1);
    vec3 disp = model_space_pos.xyz - eye_relative_group_origin;
    float distance = sqrt(dot(disp, disp));
    aabb_low = ivec3(low_x, low_y, low_z);
    aabb_high = ivec3(high_x, high_y, high_z);
    // Re-implement decide_chunk(...) == draw_raycast on GPU.
    vec3 floor_eye = floor(eye_relative_group_origin);
    vec3 aabb_nearest = clamp(floor_eye, aabb_low, aabb_high);
    disp = aabb_nearest - floor_eye;
    aabb_residue_coord = model_space_pos.xyz;
    float squared_dist = dot(disp, disp);
    bool draw_raycast = squared_dist
        == clamp(squared_dist, raycast_thresh_squared, far_plane_squared);
    // Draw as degenerate triangle if this chunk is not meant for raycasting.
    gl_Position = draw_raycast ? mvp_matrix * model_space_pos
                               : vec4(0,0,0,1);
    aabb_half_normal = unit_box_normals[gl_VertexID].xyz * 0.5;
}
