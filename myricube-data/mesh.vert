




























// Instanced inputs:
// Residue coordinates and bitfield of visible +/- x,y,z faces.
layout(location=PACKED_VERTEX_IDX) in int packed_residue_face_bits;
// Color of the voxel
layout(location=PACKED_COLOR_IDX) in int packed_color;

// Non-instanced inputs:
// Position coordinate of unit box:
layout(location=UNIT_BOX_VERTEX_IDX) in vec3 unit_box_vertex;
// Bit corresponding to the face being drawn of the unit box.  for
// example, if we're drawing the +x face, this is POS_X_FACE_BIT.
layout(location=UNIT_BOX_FACE_BIT_IDX) in int face_bit;
// Texture coordinate (for now, just used for the color border).
layout(location=UNIT_BOX_UV_IDX) in vec2 in_uv;

out vec3 v_color;
out vec3 v_residue_coord;
out vec2 v_uv;
uniform mat4 mvp_matrix;

void main() {
    // Draw as a degenerate triangle if the face we're drawing is not
    // marked visible.
    if ((face_bit & packed_residue_face_bits) == 0) {
        gl_Position = vec4(0,0,0,1);
        return;
    }

    // Reposition the unit box in the correct location.
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
    gl_Position = mvp_matrix * model_space_position;

    // Unpack the color.
    float red   = ((packed_color >> 16) & 255) * (1./255.);
    float green = ((packed_color >> 8) & 255) * (1./255.);
    float blue  = (packed_color & 255) * (1./255.);
    v_color = vec3(red, green, blue);

    v_uv = in_uv;
}
