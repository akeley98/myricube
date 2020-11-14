#version 460
#extension GL_GOOGLE_include_directive : enable

#define BORDER_WIDTH_LOW  0.075
#define BORDER_WIDTH_HIGH 0.45
#define BORDER_FADE  0.75
#define BORDER_DIST_LOW  100
#define BORDER_DIST_HIGH 350

#include "PushConstant.glsl"

layout(push_constant) uniform PushConstantBlock {
    PushConstant pc;
} pc;

layout(location=0) in vec3 v_color;
layout(location=1) in vec3 v_residue_coord;
layout(location=2) in vec2 v_uv;

layout(location=0) out vec4 out_color;

#include "fog_border.glsl"

void main() {
    vec3 disp = v_residue_coord - pc.pc.eye_relative_group_origin.xyz;
    float dist_squared = dot(disp, disp);
    vec3 fog_color = fog_color_from_world_direction(disp);
    out_color = fog_border_color(
        v_color, dist_squared, v_uv, fog_color);
}
