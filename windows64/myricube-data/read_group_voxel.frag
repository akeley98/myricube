




































// The voxel data for a given chunk group is given to the raycast
// fragment shader as a 3D texture. This shader declares the binding
// for the texture and a helper function for reading it.

uniform sampler3D chunk_group_texture;

vec4 read_group_voxel(ivec3 residue)
{
    if (clamp(residue, ivec3(0), ivec3(GROUP_SIZE - 1)) != residue) {
        return vec4(1, 1, 0, 1);
    }
    return texelFetch(chunk_group_texture, residue, 0);
}
