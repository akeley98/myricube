



















in vec3 color;
in vec3 model_space_position_;
out vec4 out_color;
void main() {
    const float d = BORDER_WIDTH;
    float x = model_space_position_.x;
    int x_border = (x - floor(x + d) < d) ? 1 : 0;
    float y = model_space_position_.y;
    int y_border = (y - floor(y + d) < d) ? 1 : 0;
    float z = model_space_position_.z;
    int z_border = (z - floor(z + d) < d) ? 1 : 0;
    float scale = (x_border + y_border + z_border >= 2) ? 0.5 : 1.0;
    out_color = vec4(color * scale, 1);
}
