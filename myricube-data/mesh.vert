





























uniform mat4 mvp_matrix;
uniform uvec3 chunk_offset_in_group;
uniform sampler3D chunk_blocks;

out vec4 v_color;
out vec3 v_residue_coord;
out vec4 mvp_x_hat;
out vec4 mvp_y_hat;
out vec4 mvp_z_hat;
out int visible_face_bits;

void main() {
    uint x_residue = uint(gl_InstanceID) % CHUNK_SIZE
                   + chunk_offset_in_group.x;
    uint y_residue = (uint(gl_InstanceID) / CHUNK_SIZE) % CHUNK_SIZE
                   + chunk_offset_in_group.y;
    uint z_residue = (uint(gl_InstanceID) / (CHUNK_SIZE * CHUNK_SIZE))
                   % CHUNK_SIZE + chunk_offset_in_group.z;

    // model space == coordinate system of the chunk group (residue coord)
    vec4 model_space_position = vec4(x_residue, y_residue, z_residue, 1);
    gl_Position = mvp_matrix * model_space_position;
    v_residue_coord = model_space_position.xyz;
    mvp_x_hat = mvp_matrix * vec4(1, 0, 0, 0);
    mvp_y_hat = mvp_matrix * vec4(0, 1, 0, 0);
    mvp_z_hat = mvp_matrix * vec4(0, 0, 1, 0);

    const float rcp = 1.0 / GROUP_SIZE;
    v_color = textureLod(chunk_blocks, model_space_position.xyz * rcp, 0);

    vec3 position_plus_x = vec3(x_residue + 1, y_residue, z_residue);
    vec3 position_minus_x = vec3(x_residue - 1, y_residue, z_residue);
    vec3 position_plus_y = vec3(x_residue, y_residue + 1, z_residue);
    vec3 position_minus_y = vec3(x_residue, y_residue - 1, z_residue);
    vec3 position_plus_z = vec3(x_residue, y_residue, z_residue + 1);
    vec3 position_minus_z = vec3(x_residue, y_residue, z_residue - 1);

    visible_face_bits = 0;
    if (v_color.a >= .5) {
        if (textureLod(chunk_blocks, position_plus_x * rcp, 0).a < .5) {
            visible_face_bits |= POS_X_FACE_BIT;
        }
        if (textureLod(chunk_blocks, position_minus_x * rcp, 0).a < .5) {
            visible_face_bits |= NEG_X_FACE_BIT;
        }
        if (textureLod(chunk_blocks, position_plus_y * rcp, 0).a < .5) {
            visible_face_bits |= POS_Y_FACE_BIT;
        }
        if (textureLod(chunk_blocks, position_minus_y * rcp, 0).a < .5) {
            visible_face_bits |= NEG_Y_FACE_BIT;
        }
        if (textureLod(chunk_blocks, position_plus_z * rcp, 0).a < .5) {
            visible_face_bits |= POS_Z_FACE_BIT;
        }
        if (textureLod(chunk_blocks, position_minus_z * rcp, 0).a < .5) {
            visible_face_bits |= NEG_Z_FACE_BIT;
        }
    }
}
