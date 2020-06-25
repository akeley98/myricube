











layout(location=PACKED_VERTEX_IDX) in int packed_vertex;
layout(location=PACKED_COLOR_IDX) in int packed_color;
out vec3 color;
out vec3 model_space_position_;
uniform mat4 mvp_matrix;
void main() {
    float x = float(packed_vertex & 255);
    float y = float((packed_vertex >> 8) & 255);
    float z = float((packed_vertex >> 16) & 255);
    vec4 model_space_position = vec4(x, y, z, 1);
    gl_Position = mvp_matrix * model_space_position;
    float red   = ((packed_color >> 16) & 255) * (1./255.);
    float green = ((packed_color >> 8) & 255) * (1./255.);
    float blue  = (packed_color & 255) * (1./255.);
    color = vec3(red, green, blue);
    model_space_position_ = model_space_position.xyz;
}
