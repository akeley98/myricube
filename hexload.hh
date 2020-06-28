// Storer/loader for a hexadecimal human-readable VoxelWorld format.
//
// The file should start with the line "myricube:hex" exactly; each
// subsequent line is for one voxel, format is
//
// [8 hex digits x],[8 hex digits y],[8 hex digits z]:[RGB, 2 hex digits each]
//
// Each line is terminated with \n (Unix/binary format)
//
// x,y,z are cast to 32-bit unsigned integers.

#ifndef MYRICUBE_HEXLOAD_HH_
#define MYRICUBE_HEXLOAD_HH_

#include "myricube.hh"

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>

#include "chunk.hh"

namespace myricube {

// Write the world to the named file. Return true if successful.
inline bool write_hex(const VoxelWorld& world, const std::string& filename)
noexcept
{
    FILE* file = fopen(filename.c_str(), "wb");
    if (file == nullptr) {
        fprintf(stderr, "Failed to open '%s'\n", filename.c_str());
        return false;
    }
    fprintf(file, "myricube:hex\n");

    bool success = true;
    int voxels_written = 0;
    world.map_voxels([file, &success, &filename, &voxels_written]
        (Voxel v, glm::ivec3 c)
    {
        assert(v.visible);
        auto x = (unsigned long)uint32_t(c.x);
        auto y = (unsigned long)uint32_t(c.y);
        auto z = (unsigned long)uint32_t(c.z);
        auto color = (unsigned long)(uint32_t(v.red) << 16
                                   | uint32_t(v.green) << 8
                                   | uint32_t(v.blue));
        int chars = fprintf(file, "%08lX,%08lX,%08lX:%06lX\n", x, y, z, color);
        if (success and chars != 34) {
            fprintf(stderr, "Failed to write all 34 chars in line to '%s' "
                "(%s)\n", filename.c_str(), strerror(errno));
            success = false;
        }
        ++voxels_written;
        if (voxels_written % (1 << 20) == 0 and success) {
            fprintf(stderr, "%iM voxels written\n", voxels_written >> 20);
        }
    });
    int err = fclose(file);
    if (err != 0 and success) {
        fprintf(stderr, "fclose for '%s' failed (%s)\n", filename.c_str(),
            strerror(errno));
        success = false;
    }
    return success;
}

// Set the voxels in the given VoxelWorld using data from the named file.
inline bool read_hex(VoxelWorld& world, const std::string& filename)
noexcept
{
    int8_t hex_digit_value_table[103] = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, -1, -1, -1, -1, -1, -1,
        -1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, 10, 11, 12, 13, 14, 15
    };

    FILE* file = fopen(filename.c_str(), "rb");
    if (file == nullptr) {
        fprintf(stderr, "Failed to open '%s' (%i) %s\n",
            filename.c_str(), errno, strerror(errno));
        return false;
    }
    bool success = true;
    const char magic[] = "myricube:hex\n";
    int line_number = 1;

    for (unsigned i = 0; i < sizeof magic - 1; ++i) {
        if (fgetc(file) != magic[i]) {
            fprintf(stderr, "Expected file '%s' to start with %s",
                filename.c_str(), magic);
            fclose(file);
            return false;
        }
    }

    auto read_hex_digits = [&] (int digits) -> int32_t
    {
        uint32_t value = 0;
        for (int i = 0; i < digits; ++i) {
            int c = fgetc(file);
            if (c == EOF and success) {
                fprintf(stderr, "Unexpected EOF %s:%i\n",
                    filename.c_str(), line_number);
                success = false;
            }
            int hex_digit = -1;
            if (unsigned(c) < sizeof hex_digit_value_table) {
                hex_digit = hex_digit_value_table[c];
            }
            if (hex_digit == -1) {
                if (success) {
                    fprintf(stderr, "%c (\\x%02x) not a hex digit %s:%i\n",
                        c, c, filename.c_str(), line_number);
                    success = false;
                }
                hex_digit = 0;
            }
            value = value << 4 | uint32_t(hex_digit);
        }
        return int32_t(value);
    };

    auto read_expected_char = [&] (int expected_char)
    {
        int c = fgetc(file);
        if (c != expected_char and success) {
            if (c == EOF) {
                fprintf(stderr, "Unexpected EOF %s:%i\n",
                    filename.c_str(), line_number);
            }
            else {
                fprintf(stderr, "Unexpected %c (\\x%02x) %s:%i\n",
                    c, c, filename.c_str(), line_number);
            }
            success = false;
        }
    };

    PositionedChunkGroup* hint = nullptr;
    int voxels_set = 0;

    while (1) {
        ++line_number;

        int c = fgetc(file);
        if (c == EOF) break;
        ungetc(c, file);

        int32_t x = read_hex_digits(8);
        read_expected_char(',');
        int32_t y = read_hex_digits(8);
        read_expected_char(',');
        int32_t z = read_hex_digits(8);
        read_expected_char(':');
        auto red = uint8_t(read_hex_digits(2));
        auto green = uint8_t(read_hex_digits(2));
        auto blue = uint8_t(read_hex_digits(2));
        read_expected_char('\n');
        
        world.set(glm::ivec3(x, y, z), Voxel(red, green, blue), &hint);
        
        ++voxels_set;
        if (voxels_set % (1 << 20) == 0 and success) {
            fprintf(stderr, "%iM voxels set\n", voxels_set >> 20);
        }
    }

    int err = fclose(file);
    if (err != 0 and success) {
        fprintf(stderr, "fclose %s failed: %s\n",
            filename.c_str(), strerror(errno));
        success = false;
    }
    return success;
}

} // end namespace

#endif /* !MYRICUBE_HEXLOAD_HH_ */
