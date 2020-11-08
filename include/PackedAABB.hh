// Just a type for computing and storing the AABB of a chunk.

#ifndef MYRICUBE_PACKEDAABB_HH_
#define MYRICUBE_PACKEDAABB_HH_

#include "myricube.hh"
#include "voxels.hh"

namespace myricube {

constexpr uint32_t x_shift = 0;
constexpr uint32_t y_shift = 8;
constexpr uint32_t z_shift = 16;

// For memory efficiency, the AABB of a chunk is stored in packed
// format on the GPU.
struct PackedAABB
{
    uint32_t packed_low;
    uint32_t packed_high;

    static_assert(group_size <= 255,
                  "group too big for 8-bit unsigned coordinates.");

    PackedAABB() = default;

    PackedAABB(glm::ivec3 aabb_low, glm::ivec3 aabb_high)
    {
        auto x = aabb_low.x;
        auto y = aabb_low.y;
        auto z = aabb_low.z;
        packed_low = uint32_t(x) << x_shift
                   | uint32_t(y) << y_shift
                   | uint32_t(z) << z_shift;
        x = aabb_high.x;
        y = aabb_high.y;
        z = aabb_high.z;
        packed_high = uint32_t(x) << x_shift
                    | uint32_t(y) << y_shift
                    | uint32_t(z) << z_shift;
    }

    // Compute the AABB (in residue coordinates) of the given chunk.
    // The residue coordinates of the lower corner of the chunk must
    // also be given (as the AABB is relative to the chunk group's
    // origin, not the chunk itself).
    PackedAABB(const BinChunk& chunk, glm::ivec3 chunk_residue)
    {
        // {xyz}_array[n] has visible_bit set iff there's some visible
        // voxel in the chunk with n == the voxel's {xyz} coordinate.
        uint32_t x_array[chunk_size] = { 0 };
        uint32_t y_array[chunk_size] = { 0 };
        uint32_t z_array[chunk_size] = { 0 };

        uint32_t all_orrd = 0;

        for (size_t z = 0; z < chunk_size; ++z) {
            for (size_t y = 0; y < chunk_size; ++y) {
                for (size_t x = 0; x < chunk_size; ++x) {
                    uint32_t voxel = chunk.voxel_array[z][y][x];
                    all_orrd |= voxel;
                    x_array[x] |= voxel;
                    y_array[y] |= voxel;
                    z_array[z] |= voxel;
                }
            }
        }

        // Return trivial PackedAABB if chunk is empty.
        if ((all_orrd & visible_bit) == 0) {
            packed_low = 0;
            packed_high = 0;
            return;
        }

        // Otherwise, we actually have to calculate the AABB. Do it
        // in the coordinate system of the chunk first.
        glm::ivec3 chunk_aabb_low(chunk_size, chunk_size, chunk_size);
        glm::ivec3 chunk_aabb_high(0, 0, 0);

        for (int32_t i = 0; i < chunk_size; ++i) {
            if (x_array[i] & visible_bit) {
                chunk_aabb_low.x = glm::min(chunk_aabb_low.x, i);
                chunk_aabb_high.x = glm::max(chunk_aabb_low.x, i+1);
            }
            if (y_array[i] & visible_bit) {
                chunk_aabb_low.y = glm::min(chunk_aabb_low.y, i);
                chunk_aabb_high.y = glm::max(chunk_aabb_low.y, i+1);
            }
            if (z_array[i] & visible_bit) {
                chunk_aabb_low.z = glm::min(chunk_aabb_low.z, i);
                chunk_aabb_high.z = glm::max(chunk_aabb_low.z, i+1);
            }
        }

        // Reposition computed AABB by the chunk's residue coordinate.
        *this = PackedAABB(
            chunk_aabb_low + chunk_residue,
            chunk_aabb_high + chunk_residue);
    }

    uint8_t low_x() const
    {
        return uint8_t((packed_low >> x_shift) & 255);
    }
    uint8_t low_y() const
    {
        return uint8_t((packed_low >> y_shift) & 255);
    }
    uint8_t low_z() const
    {
        return uint8_t((packed_low >> z_shift) & 255);
    }
    uint8_t high_x() const
    {
        return uint8_t((packed_high >> x_shift) & 255);
    }
    uint8_t high_y() const
    {
        return uint8_t((packed_high >> y_shift) & 255);
    }
    uint8_t high_z() const
    {
        return uint8_t((packed_high >> z_shift) & 255);
    }
};

} // end namespace myricube

#endif /* !MYRICUBE_PackedAABB_HH_ */
