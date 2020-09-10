




































// The voxel data for a given chunk group is given to the raycast
// fragment shader as a GROUP_SIZE x GROUP_SIZE x GROUP_SIZE 3D array
// (in [z][y][x] order) of packed 32-bit ints (defines visibility,
// red, green, blue), as well as a packed bit array that only carries
// visibility information (this helps speed up the common case of
// raycasting through empty space).
//
// This shader defines the SSBO and UBO bindings that carries this
// info, and a helper function for unpacking it.

// 3D array of voxels, in [z][y][x] order.
layout(std430, binding=CHUNK_GROUP_VOXELS_PROGRAM_INDEX)
readonly buffer ChunkGroupVoxels {
    uint voxel_colors[GROUP_SIZE][GROUP_SIZE][GROUP_SIZE];
} chunk_group_voxels;

// Enumerate voxels in a chunk group in z, y, x order (i.e. x varies
// fastest), then voxel number i is visible iff bit i is nonzero.
// Bit 0 is the 1s place of bits[0].
layout(packed, binding=CHUNK_GROUP_VISIBLE_BITS_PROGRAM_INDEX)
uniform ChunkGroupVisibleBits {
    uint bits[GROUP_SIZE * GROUP_SIZE * GROUP_SIZE / 32];
} chunk_group_visible;
// I have to use packed instead of std140 to avoid the 16-byte stride
// behavior. Ballsey but what else can I do???

// If true we skip checking the above UBO.
uniform bool disable_chunk_group_visible_bits;

vec4 read_group_voxel(ivec3 residue)
{
    if (clamp(residue, ivec3(0), ivec3(GROUP_SIZE - 1)) != residue) {
        return vec4(1, 1, 0, 1);
    }

    // First check the visible bits, and exit early if it's not visible.
    if (!disable_chunk_group_visible_bits) {
        int bit_index = residue.z;
        bit_index = residue.y + GROUP_SIZE * bit_index;
        bit_index = residue.x + GROUP_SIZE * bit_index;
        int array_index = bit_index >> 5;
        bit_index &= 31;
        uint visible =
            1u & (chunk_group_visible.bits[array_index] >> bit_index);
        if (visible == 0) return vec4(0);
    }

    uint packed_color = chunk_group_voxels.voxel_colors
        [residue.z][residue.y][residue.x];

    vec4 color;
    color.a = (packed_color & uint(VISIBLE_BIT)) != 0 ? 1.0 : 0.0;
    color.r = ((packed_color >> RED_SHIFT) & 255u) * (1 / 255.);
    color.g = ((packed_color >> GREEN_SHIFT) & 255u) * (1 / 255.);
    color.b = ((packed_color >> BLUE_SHIFT) & 255u) * (1 / 255.);
    if (disable_chunk_group_visible_bits) return color;

    // For now, check if the bitfield is correct.
    return color.a == 0 ? vec4(1) : color;
}
