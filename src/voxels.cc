// Implementation functions for WorldHandle.  Basically we just have
// to negotiate with the operating system to get memory-mapped views
// of chunk groups on disk (or of the in-memory filesystem).

#include "voxels.hh"

#include <errno.h>
#include <mutex>
#include <new>
#include <stdexcept>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

#include "EnvVar.hh"
#include "map.hh"

namespace myricube {

void unmap_mut_bin_chunk_group(BinChunkGroup* ptr)
{
    unmap(ptr);
}

void unmap_bin_chunk_group(const BinChunkGroup* ptr)
{
    unmap(ptr);
}

void unmap_bin_group_bitfield(BinGroupBitfield* ptr)
{
    unmap(ptr);
}

WorldHandle::WorldHandle(const filename_string& arg)
{
    auto on_bad_filename = [&arg]
    {
#ifdef MYRICUBE_WINDOWS
        fprintf(stderr, "Got bad filename: '%ls'\n", arg.c_str());
#endif
        throw std::runtime_error("Expected world file name to be "
            + std::string(world_filename) + " or end in /"
#ifndef MYRICUBE_WINDOWS
            + "\nGot " + arg
#endif
        );
    };

    // First we need to recover the directory name for the world.
    // Remove until we get to a trailing forward slash, or, on
    // Windows only, a backwards slash.
    directory_trailing_slash = arg;
    while (1) {
        if (directory_trailing_slash.empty()) {
            on_bad_filename();
        }
        if (directory_trailing_slash.back() == filename_char('/')) break;
#ifdef MYRICUBE_WINDOWS
        if (directory_trailing_slash.back() == filename_char('\\')) break;
#endif
        directory_trailing_slash.pop_back();
    }

    auto bin_world_filename = filename_concat_c_str(
        directory_trailing_slash, world_filename);

    // Check that the file name was as expected (world.myricube at the moment),
    // or the argument was a directory with a trailing slash.
    if (bin_world_filename != arg and arg != directory_trailing_slash) {
        on_bad_filename();
    }

    // Now I need to memory-map the BinWorld file, and get the shared_ptr
    // to unmap it when finished.
    bin_world = std::shared_ptr<BinWorld>(
        map_file<BinWorld>(bin_world_filename, create_flag),
        [](BinWorld* p_world) { unmap(p_world); });

    // Finally check the magic number.
    // Ignore the endian bit as that will be checked later.
    auto xord = bin_world->magic_number ^ bin_world->expected_magic;
    if ((xord & ~new_endian_magic) != 0) {
        #ifdef MYRICUBE_WINDOWS
            fwprintf(stderr,
                L"Incorrect magic number: %ls",
                bin_world_filename.c_str());
            throw std::runtime_error("Incorrect magic number");
        #else
            throw std::runtime_error("Incorrect magic number: "
                + bin_world_filename);
        #endif
    }
}

// Filename (not including directory) of the file storing the named
// chunk group.
std::string group_coord_filename(glm::ivec3 group_coord)
{
    std::string result;
    result.resize(50);
    int bytes = sprintf(
        result.data(),
        "z%08X-%08X-%08X",
        unsigned(group_coord.x),
        unsigned(group_coord.y),
        unsigned(group_coord.z));
    assert(bytes < int(result.size()));
    result.resize(bytes);
    return result;
}

// Older versions of myricube used a different endianness for files.
// I detected this by changing the magic number for new endianness
// files.  I detect the old magic number here and fix the endianness
// and magic number in-place if needed, but only if the user opted in
// with a nonzero myricube_endian_fix environment variable (as this is
// a potentially dangerous operation if interrupted).
//
// Return value: Return true if the endianness was wrong and we
// corrected it.
//
// NOTE: This function is not actually const of course but it doesn't
// matter since the actual files should all be read/write (they have
// to be due to the dirty bits) and this is a hack anyway.
static inline bool maybe_fix_endian(
    const filename_string& filename,
    const BinChunkGroup& group_const)
{
    BinChunkGroup& group = const_cast<BinChunkGroup&>(group_const);

    auto old_magic_number = group.expected_magic & ~new_endian_magic;
    if (group.magic_number != old_magic_number) {
        return false;
    }

    thread_local EnvVar64 endian_fix_enabled("myricube_endian_fix", 0);

    if (!endian_fix_enabled) {
        #ifdef MYRICUBE_WINDOWS
            fprintf(stderr, "Incorrect endianness: %ls\n", filename.c_str());
            throw std::runtime_error("File had incorrect endianness"
                "\nrun with environment variable myricube_endian_fix=1"
                "\nto fix in-place (backup first!)");
        #else
            throw std::runtime_error("Incorrect endianness: " + filename
                + "\nrun with environment variable myricube_endian_fix=1"
                  "\nto fix in-place (backup first!)");
        #endif
    }

    // I just realized the horrible things this function can do if run
    // concurrently, so fix this here. This won't protect in case
    // multiple myricube instances are running, but it's better than
    // nothing (this function truly is a HORRIBLE hack, but I kinda
    // need it now for bug-for-bug compatibility).
    static std::mutex the_mutex;
    std::lock_guard guard(the_mutex);
    if (group.magic_number != old_magic_number) {
        fprintf(stderr, "maybe_fix_endian: fixed by another thread\n");
        return false;
    }

    for (int zH = 0; zH < edge_chunks; ++zH) {
    for (int yH = 0; yH < edge_chunks; ++yH) {
    for (int xH = 0; xH < edge_chunks; ++xH) {
        BinChunk& chunk = group.chunk_array[zH][yH][xH];

        for (int zL = 0; zL < chunk_size; ++zL) {
        for (int yL = 0; yL < chunk_size; ++yL) {
        for (int xL = 0; xL < chunk_size; ++xL) {
            uint32_t* p_voxel = &chunk.voxel_array[zL][yL][xL];
            uint32_t old_voxel = *p_voxel;
            uint32_t new_voxel = 0;
            new_voxel |= ((old_voxel >> 0) & 0xFF) << 24;
            new_voxel |= ((old_voxel >> 8) & 0xFF) << 16;
            new_voxel |= ((old_voxel >> 16) & 0xFF) << 8;
            new_voxel |= ((old_voxel >> 24) & 0xFF) << 0;
            *p_voxel = new_voxel;
        }
        }
        }
    }
    }
    }

    group.magic_number |= new_endian_magic;
    #ifdef MYRICUBE_WINDOWS
        fprintf(stderr, "Fixed endianness of %ls\n", filename.c_str());
    #else
        fprintf(stderr, "Fixed endianness of %s\n", filename.c_str());
    #endif
    return true;
}

UPtrMutChunkGroup WorldHandle::mut_chunk_group(glm::ivec3 group_coord)
{
    // Map the correct file for this chunk group.
    auto filename = filename_concat_c_str(
        directory_trailing_slash, group_coord_filename(group_coord).c_str());
    int flags = create_flag;
    auto result =
        UPtrMutChunkGroup(map_file<BinChunkGroup>(filename, &flags));
    assert(result != nullptr);

    // Unconditionally set the bit in the BinGroupBitfield for this
    // chunk group, even if it was already created. This is more
    // robust against data races and I expect WorldMutator to avoid
    // calling us as much as possible.
    bitfield_for_chunk_group(group_coord)->set_chunk_group_on_disk(group_coord);

    // Check magic number and return.
    if (result->magic_number != result->expected_magic) {
        bool okay = maybe_fix_endian(filename, *result);
        if (okay) goto its_okay;
        #ifdef MYRICUBE_WINDOWS
            fwprintf(stderr,
                L"Incorrect magic number: %ls",
                filename.c_str());
            throw std::runtime_error("Incorrect magic number");
        #else
            throw std::runtime_error("Incorrect magic number: " + filename);
        #endif
    }
  its_okay:
    return result;
}

UPtrChunkGroup WorldHandle::view_chunk_group(glm::ivec3 group_coord) const
{
    // Note: we don't check the group bitfield; someone else could
    // (for performance) have already checked that before calling us.
    auto filename = filename_concat_c_str(
        directory_trailing_slash, group_coord_filename(group_coord).c_str());
    auto result = UPtrChunkGroup(
        map_file<const BinChunkGroup>(filename, readonly_flag));
    if (result == nullptr) return nullptr;

    if (result->magic_number != result->expected_magic) {
        bool okay = maybe_fix_endian(filename, *result);
        if (okay) goto its_okay;
        #ifdef MYRICUBE_WINDOWS
            fwprintf(stderr,
                L"Incorrect magic number: %ls\n", filename.c_str());
            throw std::runtime_error("Incorrect magic number");
        #else
            throw std::runtime_error("Incorrect magic number: " + filename);
        #endif
    }
  its_okay:
    return result;
}

UPtrGroupBitfield
WorldHandle::bitfield_for_chunk_group(glm::ivec3 group_coord) const
{
    // Make the filename for this bitfield by appending a 'b' to
    // the filename for the (canonical) chunk group within.
    auto canonical = BinGroupBitfield::canonical_group_coord(group_coord);
    std::string n = group_coord_filename(canonical) + "b";
    auto filename = filename_concat_c_str(
        directory_trailing_slash, n.c_str());

    auto result =
        UPtrGroupBitfield(map_file<BinGroupBitfield>(filename, create_flag));
    assert(result != nullptr);
    return result;
}

} // end namespace myricube
