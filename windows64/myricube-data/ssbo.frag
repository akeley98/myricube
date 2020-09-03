




































// The voxel data for a given chunk group is given to the raycast
// fragment shader as a GROUP_SIZE x GROUP_SIZE x GROUP_SIZE 3D array
// (in [z][y][x] order) of packed 32-bit ints (defines visibility,
// red, green, blue). This shader defines the SSBO binding that
// carries this info, and a helper function for unpacking it.

layout(std430, binding=CHUNK_GROUP_VOXELS_PROGRAM_INDEX)
buffer ChunkGroupVoxels {
    int voxel_colors[GROUP_SIZE][GROUP_SIZE][GROUP_SIZE];
} chunk_group_voxels;

vec4 read_group_voxel(ivec3 residue)
{
    if (clamp(residue, ivec3(0), ivec3(GROUP_SIZE - 1)) != residue) {
        return vec4(1, 1, 0, 1);
    }
    int packed_color = chunk_group_voxels.voxel_colors
        [residue.z][residue.y][residue.x];

    vec4 color;
    color.a = (packed_color & VISIBLE_BIT) != 0 ? 1.0 : 0.0;
    color.r = ((packed_color >> RED_SHIFT) & 255) * (1 / 255.);
    color.g = ((packed_color >> GREEN_SHIFT) & 255) * (1 / 255.);
    color.b = ((packed_color >> BLUE_SHIFT) & 255) * (1 / 255.);
    return color;
}
