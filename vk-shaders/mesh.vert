#version 460
#extension GL_GOOGLE_include_directive : require

#include "meshcommon.glsl"
#include "PushConstant.glsl"
#include "srgb.glsl"

// Instanced inputs:
// Residue coordinates and bitfield of visible +/- x,y,z faces.
layout(location=0) in uint packed_residue_face_bits;
// Color of the voxel
layout(location=1) in uint packed_color;

layout(location=0) out vec3 v_color;
layout(location=1) out vec3 v_residue_coord;
layout(location=2) out vec2 v_uv;

// Non-instanced inputs: unit box vertices built-into vertex shader.
// Draw as GL_TRIANGLES with 36 vertices. (I can't use a triangle
// strip due to the face visibility bits).

void main() {
    int face_bit = unit_box_face_bits[gl_VertexIndex];

    // Draw as a degenerate triangle if the face we're drawing is not
    // marked visible.
    if ((face_bit & packed_residue_face_bits) == 0) {
        gl_Position = vec4(0,0,0,1);
        return;
    }

    // Reposition the unit box in the correct location.
    vec3 unit_box_vertex = unit_box_verts[gl_VertexIndex];
    float x = float((packed_residue_face_bits >> X_SHIFT) & 255u)
            + unit_box_vertex.x;
    float y = float((packed_residue_face_bits >> Y_SHIFT) & 255u)
            + unit_box_vertex.y;
    float z = float((packed_residue_face_bits >> Z_SHIFT) & 255u)
            + unit_box_vertex.z;
    vec4 model_space_position = vec4(x, y, z, 1);
    v_residue_coord = model_space_position.xyz;

    // Perspective transformation. Note that the location of the chunk
    // group we're in is taken care of by the 'm' in mvp.
    gl_Position = pc.pc.mvp * model_space_position;

    // Unpack the color.
    float red   = srgb_from_u8((packed_color >> RED_SHIFT) & 255);
    float green = srgb_from_u8((packed_color >> GREEN_SHIFT) & 255);
    float blue  = srgb_from_u8((packed_color >> BLUE_SHIFT) & 255);
    v_color = vec3(red, green, blue);

    v_uv = uv_array[gl_VertexIndex];
}
