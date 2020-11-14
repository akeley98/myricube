#version 460
#extension GL_GOOGLE_include_directive : require

#define POS_X_FACE_BIT (1 << 24)
#define NEG_X_FACE_BIT (1 << 25)
#define POS_Y_FACE_BIT (1 << 26)
#define NEG_Y_FACE_BIT (1 << 27)
#define POS_Z_FACE_BIT (1 << 28)
#define NEG_Z_FACE_BIT (1 << 29)

#define X_SHIFT 0
#define Y_SHIFT 8
#define Z_SHIFT 16

#define RED_SHIFT 24
#define GREEN_SHIFT 16
#define BLUE_SHIFT 8

#include "PushConstant.glsl"

layout (push_constant) uniform PushConstantBlock {
    PushConstant pc;
} pc;

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

// Position coordinate of unit box:
vec3 unit_box_verts[36] = vec3[] (
    vec3(0, 1, 1),    vec3(0, 1, 0),    vec3(0, 0, 1),
    vec3(0, 1, 0),    vec3(0, 0, 0),    vec3(0, 0, 1), // -x face

    vec3(1, 0, 0),    vec3(1, 1, 0),    vec3(1, 0, 1),
    vec3(1, 1, 0),    vec3(1, 1, 1),    vec3(1, 0, 1), // +x face

    vec3(0, 0, 0),    vec3(1, 0, 1),    vec3(0, 0, 1),
    vec3(0, 0, 0),    vec3(1, 0, 0),    vec3(1, 0, 1), // -y face

    vec3(1, 1, 1),    vec3(1, 1, 0),    vec3(0, 1, 0),
    vec3(0, 1, 1),    vec3(1, 1, 1),    vec3(0, 1, 0), // +y face

    vec3(0, 0, 0),    vec3(0, 1, 0),    vec3(1, 0, 0),
    vec3(1, 1, 0),    vec3(1, 0, 0),    vec3(0, 1, 0), // -z face

    vec3(0, 1, 1),    vec3(1, 0, 1),    vec3(1, 1, 1),
    vec3(1, 0, 1),    vec3(0, 1, 1),    vec3(0, 0, 1));// +z face

// Bit corresponding to the face being drawn of the unit box; for
// example, if we're drawing the +x face, this is POS_X_FACE_BIT.
int unit_box_face_bits[36] = int[] (
    NEG_X_FACE_BIT,    NEG_X_FACE_BIT,    NEG_X_FACE_BIT,
    NEG_X_FACE_BIT,    NEG_X_FACE_BIT,    NEG_X_FACE_BIT,

    POS_X_FACE_BIT,    POS_X_FACE_BIT,    POS_X_FACE_BIT,
    POS_X_FACE_BIT,    POS_X_FACE_BIT,    POS_X_FACE_BIT,

    NEG_Y_FACE_BIT,    NEG_Y_FACE_BIT,    NEG_Y_FACE_BIT,
    NEG_Y_FACE_BIT,    NEG_Y_FACE_BIT,    NEG_Y_FACE_BIT,

    POS_Y_FACE_BIT,    POS_Y_FACE_BIT,    POS_Y_FACE_BIT,
    POS_Y_FACE_BIT,    POS_Y_FACE_BIT,    POS_Y_FACE_BIT,

    NEG_Z_FACE_BIT,    NEG_Z_FACE_BIT,    NEG_Z_FACE_BIT,
    NEG_Z_FACE_BIT,    NEG_Z_FACE_BIT,    NEG_Z_FACE_BIT,

    POS_Z_FACE_BIT,    POS_Z_FACE_BIT,    POS_Z_FACE_BIT,
    POS_Z_FACE_BIT,    POS_Z_FACE_BIT,    POS_Z_FACE_BIT);

// Texture coordinate (for now, just used for the color border).
vec2 uv_array[36] = vec2[] (
    vec2(1, 1),    vec2(1, 0),    vec2(0, 1),
    vec2(1, 0),    vec2(0, 0),    vec2(0, 1), // -x face

    vec2(0, 0),    vec2(1, 0),    vec2(0, 1),
    vec2(1, 0),    vec2(1, 1),    vec2(0, 1), // +x face

    vec2(0, 0),    vec2(1, 1),    vec2(0, 1),
    vec2(0, 0),    vec2(1, 0),    vec2(1, 1), // -y face

    vec2(1, 1),    vec2(1, 0),    vec2(0, 0),
    vec2(0, 1),    vec2(1, 1),    vec2(0, 0), // +y face

    vec2(0, 0),    vec2(0, 1),    vec2(1, 0),
    vec2(1, 1),    vec2(1, 0),    vec2(0, 1), // -z face

    vec2(0, 1),    vec2(1, 0),    vec2(1, 1),
    vec2(1, 0),    vec2(0, 1),    vec2(0, 0));// +z face

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
    float x = float((packed_residue_face_bits >> X_SHIFT) & 255)
            + unit_box_vertex.x;
    float y = float((packed_residue_face_bits >> Y_SHIFT) & 255)
            + unit_box_vertex.y;
    float z = float((packed_residue_face_bits >> Z_SHIFT) & 255)
            + unit_box_vertex.z;
    vec4 model_space_position = vec4(x, y, z, 1);
    v_residue_coord = model_space_position.xyz;

    // Perspective transformation. Note that the location of the chunk
    // group we're in is taken care of by the 'm' in mvp.
    gl_Position = pc.pc.mvp * model_space_position;

    // Unpack the color.
    float red   = ((packed_color >> RED_SHIFT) & 255) * (1./255.);
    float green = ((packed_color >> GREEN_SHIFT) & 255) * (1./255.);
    float blue  = ((packed_color >> BLUE_SHIFT) & 255) * (1./255.);
    v_color = vec3(red, green, blue);

    v_uv = uv_array[gl_VertexIndex];
}
