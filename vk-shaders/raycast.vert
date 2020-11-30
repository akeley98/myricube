#version 460
#extension GL_GOOGLE_include_directive : require

// Vertex shader that expands packed AABB (defines lower and upper
// corner of AABB -- passed as instanced vertex attribute) into actual
// AABBs by stretching and positioning a built-in unit box.  Draw as a
// 14-vertex triangle strip, with FIRST vertex as provoking vertex and
// counter-clockwise front (unculled) faces. However, boxes that are
// close enough to have been drawn by the mesh shader are discarded.

// https://stackoverflow.com/questions/28375338/cube-using-single-gl-triangle-strip

#define X_SHIFT 0
#define Y_SHIFT 8
#define Z_SHIFT 16

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

// Corners of the AABB to draw, in residue coordinates, packed in 8/8/8 bits.
layout(location=0) in int packed_aabb_low;
layout(location=1) in int packed_aabb_high;

#include "PushConstant.glsl"

layout(location=0)      out vec3 aabb_residue_coord;
layout(location=1) flat out ivec3 aabb_low;
layout(location=2) flat out ivec3 aabb_high;
layout(location=3) flat out vec3 aabb_half_normal;

void main() {
    int low_x = (packed_aabb_low >> X_SHIFT) & 255;
    int low_y = (packed_aabb_low >> Y_SHIFT) & 255;
    int low_z = (packed_aabb_low >> Z_SHIFT) & 255;
    int high_x = (packed_aabb_high >> X_SHIFT) & 255;
    int high_y = (packed_aabb_high >> Y_SHIFT) & 255;
    int high_z = (packed_aabb_high >> Z_SHIFT) & 255;
    vec3 f_aabb_low = vec3(low_x, low_y, low_z);
    vec3 sz = vec3(high_x, high_y, high_z) - f_aabb_low;
    vec3 unit_box_vertex = unit_box_verts[gl_VertexIndex].xyz;
    vec4 model_space_pos = vec4(unit_box_vertex * sz + f_aabb_low, 1);
    vec3 disp = model_space_pos.xyz - pc.pc.eye_relative_group_origin.xyz;
    float distance = sqrt(dot(disp, disp));
    aabb_low = ivec3(low_x, low_y, low_z);
    aabb_high = ivec3(high_x, high_y, high_z);
    // Re-implement decide_chunk(...) == draw_raycast on GPU.
    vec3 floor_eye = floor(pc.pc.eye_relative_group_origin.xyz);
    vec3 aabb_nearest = clamp(floor_eye, aabb_low, aabb_high);
    disp = aabb_nearest - floor_eye;
    aabb_residue_coord = model_space_pos.xyz;
    float squared_dist = dot(disp, disp);
    bool draw_raycast = squared_dist
        == clamp(squared_dist,
                 pc.pc.raycast_thresh_squared,
                 pc.pc.far_plane_squared);
    // Draw as degenerate triangle if this chunk is not meant for raycasting.
    gl_Position = draw_raycast ? pc.pc.mvp * model_space_pos
                               : vec4(0,0,0,1);
    aabb_half_normal = unit_box_normals[gl_VertexIndex].xyz * 0.5;
}
