





























layout (points) in;
layout (triangle_strip, max_vertices = 24) out;

in vec4 v_color[];
in vec3 v_residue_coord[];
in vec4 mvp_x_hat[];
in vec4 mvp_y_hat[];
in vec4 mvp_z_hat[];
in int visible_face_bits[];

out vec3 color;
out vec3 residue_coord;
out vec2 uv;

void emit_v000() {
    vec4 mvp_position = gl_in[0].gl_Position;
    gl_Position = mvp_position;
    color = v_color[0].rgb;
    residue_coord = v_residue_coord[0];
    EmitVertex();
}

void emit_v001() {
    vec4 mvp_position = gl_in[0].gl_Position;
    gl_Position = mvp_position + mvp_z_hat[0];
    color = v_color[0].rgb;
    residue_coord = v_residue_coord[0] + vec3(0, 0, 1);
    EmitVertex();
}

void emit_v010() {
    vec4 mvp_position = gl_in[0].gl_Position;
    gl_Position = mvp_position + mvp_y_hat[0];
    color = v_color[0].rgb;
    residue_coord = v_residue_coord[0] + vec3(0, 1, 0);
    EmitVertex();
}

void emit_v011() {
    vec4 mvp_position = gl_in[0].gl_Position;
    gl_Position = mvp_position + mvp_y_hat[0] + mvp_z_hat[0];
    color = v_color[0].rgb;
    residue_coord = v_residue_coord[0] + vec3(0, 1, 1);
    EmitVertex();
}

void emit_v100() {
    vec4 mvp_position = gl_in[0].gl_Position;
    gl_Position = mvp_position + mvp_x_hat[0];
    color = v_color[0].rgb;
    residue_coord = v_residue_coord[0] + vec3(1, 0, 0);
    EmitVertex();
}

void emit_v101() {
    vec4 mvp_position = gl_in[0].gl_Position;
    gl_Position = mvp_position + mvp_x_hat[0] + mvp_z_hat[0];
    color = v_color[0].rgb;
    residue_coord = v_residue_coord[0] + vec3(1, 0, 1);
    EmitVertex();
}

void emit_v110() {
    vec4 mvp_position = gl_in[0].gl_Position;
    gl_Position = mvp_position + mvp_x_hat[0] + mvp_y_hat[0];
    color = v_color[0].rgb;
    residue_coord = v_residue_coord[0] + vec3(1, 1, 0);
    EmitVertex();
}

void emit_v111() {
    vec4 mvp_position = gl_in[0].gl_Position;
    gl_Position = mvp_position + mvp_x_hat[0] + mvp_y_hat[0] + mvp_z_hat[0];
    color = v_color[0].rgb;
    residue_coord = v_residue_coord[0] + vec3(1, 1, 1);
    EmitVertex();
}

void main() {
    vec4 mvp_position = gl_in[0].gl_Position;

    if ((visible_face_bits[0] & POS_X_FACE_BIT) != 0) {
        uv = vec2(0, 0);
        emit_v100();
        uv = vec2(1, 0);
        emit_v110();
        uv = vec2(0, 1);
        emit_v101();
        uv = vec2(1, 1);
        emit_v111();
        EndPrimitive();
    }

    if ((visible_face_bits[0] & NEG_X_FACE_BIT) != 0) {
        uv = vec2(0, 0);
        emit_v000();
        uv = vec2(1, 0);
        emit_v001();
        uv = vec2(0, 1);
        emit_v010();
        uv = vec2(1, 1);
        emit_v011();
        EndPrimitive();
    }

    if ((visible_face_bits[0] & POS_Y_FACE_BIT) != 0) {
        uv = vec2(0, 0);
        emit_v010();
        uv = vec2(1, 0);
        emit_v011();
        uv = vec2(0, 1);
        emit_v110();
        uv = vec2(1, 1);
        emit_v111();
        EndPrimitive();
    }

    if ((visible_face_bits[0] & NEG_Y_FACE_BIT) != 0) {
        uv = vec2(0, 0);
        emit_v000();
        uv = vec2(1, 0);
        emit_v100();
        uv = vec2(0, 1);
        emit_v001();
        uv = vec2(1, 1);
        emit_v101();
        EndPrimitive();
    }

    if ((visible_face_bits[0] & POS_Z_FACE_BIT) != 0) {
        uv = vec2(0, 0);
        emit_v001();
        uv = vec2(1, 0);
        emit_v101();
        uv = vec2(0, 1);
        emit_v011();
        uv = vec2(1, 1);
        emit_v111();
        EndPrimitive();
    }

    if ((visible_face_bits[0] & NEG_Z_FACE_BIT) != 0) {
        uv = vec2(0, 0);
        emit_v000();
        uv = vec2(1, 0);
        emit_v010();
        uv = vec2(0, 1);
        emit_v100();
        uv = vec2(1, 1);
        emit_v110();
        EndPrimitive();
    }
}
