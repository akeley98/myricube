















layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

in vec3 v_color[];
in vec3 v_residue_coord[];

out vec3 color;
out vec3 residue_coord;

void main() {
    gl_Position = gl_in[0].gl_Position;
    color = v_color[0];
    residue_coord = v_residue_coord[0];
    EmitVertex();

    gl_Position = gl_in[1].gl_Position;
    // color = v_color[1].rgb * 0.5;
    color = v_color[1].rgb;
    residue_coord = v_residue_coord[1];
    EmitVertex();

    gl_Position = gl_in[2].gl_Position;
    color = v_color[2];
    residue_coord = v_residue_coord[2];
    EmitVertex();

    EndPrimitive();
}
