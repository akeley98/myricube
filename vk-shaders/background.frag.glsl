#version 460
#extension GL_GOOGLE_include_directive : enable

layout(location=0) in vec2 normalized_screen_xy;
layout(location=0) out vec4 color;

#include "PushConstant.glsl"
#include "fog_border.glsl"

// Basically, reverse engineer this fragment's position and compare it
// to the eye position to get a direction-we're-looking at, which can
// be converted to a fog color.
//
// HACK to reuse the same push constant layout that the fog functions
// expect, MVP really means inverse view-projection matrix, and
// eye_relative_group_origin is the eye's residue coord in the chunk
// group it's currently in.
void main() {
    mat4 inverse_vp = pc.pc.mvp;
    vec4 v = inverse_vp * vec4(normalized_screen_xy, 1.0, 1.0);
    vec3 frag_world_position = v.xyz / v.w;
    vec3 direction = frag_world_position - pc.pc.eye_relative_group_origin.xyz;
    vec3 fog_color = fog_color_from_world_direction(direction);
    color = vec4(fog_color, 1.0);
}
