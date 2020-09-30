// Implementation functions for WorldHandle.
// Basically we just have to negotiate with the operating system to
// get memory-mapped views of chunk groups on disk.

#include "voxels.hh"

#include <errno.h>
#include <new>
#include <stdexcept>
#include <string.h>
#include <type_traits>

#ifdef MYRICUBE_WINDOWS
// todo
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

namespace myricube {

// Given the name of a file, mmap it and return the pointer. If the
// file doesn't exist or has 0 size, either return nullptr or create
// it depending on flags (defined below). If created, the T
// constructor is run in-place on the file. THIS IS EXPECTED TO NOT
// THROW AN EXCEPTION.
//
// Throw an exception on failure to open for any reason other than
// file not existing. An exception is also thrown if the file size is
// nonzero and doesn't match that of type T (for now I don't
// anticipate mapping arrays).
template <typename T> T* map_file(const filename_string& filename, int* flags);
template <typename T> T* map_file(const filename_string& filename, int flags)
{
    return map_file<T>(filename, &flags);
}

// Forget about readonly_flag for now because of mutable atomic counter.
// constexpr int readonly_flag = 1; // Added automatically if const ptr requested.

constexpr int create_flag = 2;   // If set, create file if needed.
constexpr int file_created_flag = 4; // Set by the map_file function iff
                                     // the file needed to be created.

void* map_file_impl(const filename_string&, size_t sz, int* flags);

template <typename T> T* map_file(const filename_string& filename, int* flags)
{
    // if (std::is_const_v<T>) *flags |= readonly_flag;
    auto ptr = static_cast<T*>(map_file_impl(filename, sizeof(T), flags));
    if (*flags & file_created_flag) {
        new((void*)ptr) T;
    }
    // TODO think about potential race / interruption condition:
    // create file in separate location and move-in when ready?
    return ptr;
}

template <typename T> void unmap_file(const T* ptr)
{
    // fprintf(stderr, "unmap_file of size %ld\n", long(sizeof(T)));
    auto code = munmap((void*)ptr, sizeof(T));
    assert(code == 0); // Seems like munmap only fails due to logic errors.
}

void unmap_mut_bin_chunk_group(BinChunkGroup* ptr)
{
    unmap_file(ptr);
}

void unmap_bin_chunk_group(const BinChunkGroup* ptr)
{
    unmap_file(ptr);
}

WorldHandle::WorldHandle(const filename_string& world_filename_arg)
{
    auto on_bad_filename = [&world_filename_arg]
    {
        throw std::runtime_error("Expected world file name to be "
            + std::string(world_filename)
#ifndef MYRICUBE_WINDOWS
            + "\nGot " + world_filename_arg
#endif
        );
    };

    // First we need to recover the directory name for the world.
    // Remove until we get to a trailing forward slash, or, on
    // Windows only, a backwards slash.
    directory_trailing_slash = world_filename_arg;
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

    // Check that the file name was as expected (world.myricube at the moment).
    filename_string check_string =
        filename_concat_c_str(directory_trailing_slash, world_filename);
    if (check_string != world_filename_arg) on_bad_filename();

    // Now I need to memory-map the BinWorld file, and get the shared_ptr
    // to unmap it when finished.
    auto bin_world_filename = filename_concat_c_str(
        directory_trailing_slash, world_filename);
    bin_world = std::shared_ptr<BinWorld>(
        map_file<BinWorld>(bin_world_filename, create_flag),
        unmap_file<BinWorld>);

    // Finally check the magic number.
    if (bin_world->magic_number != bin_world->expected_magic) {
        throw std::runtime_error("Incorrect magic number: "
            + bin_world_filename);
    }
}

std::string group_coord_filename(glm::ivec3 group_coord)
{
    std::string result;
    result.resize(50);
    int bytes = sprintf(
        result.data(),
        "%08X-%08X-%08X",
        unsigned(group_coord.x),
        unsigned(group_coord.y),
        unsigned(group_coord.z));
    assert(bytes < int(result.size()));
    return result;
}

UPtrMutChunkGroup WorldHandle::mut_chunk_group(glm::ivec3 group_coord)
{
    auto filename = filename_concat_c_str(
        directory_trailing_slash, group_coord_filename(group_coord).c_str());
    auto result =
        UPtrMutChunkGroup(map_file<BinChunkGroup>(filename, create_flag));
    assert(result != nullptr);
    if (result->magic_number != result->expected_magic) {
        throw std::runtime_error("Incorrect magic number: " + filename);
    }
    return result;
}

UPtrChunkGroup WorldHandle::view_chunk_group(glm::ivec3 group_coord) const
{
    auto filename = filename_concat_c_str(
        directory_trailing_slash, group_coord_filename(group_coord).c_str());
    auto result = UPtrChunkGroup(map_file<const BinChunkGroup>(filename, 0));
    if (result == nullptr) return nullptr;
    if (result->magic_number != result->expected_magic) {
        throw std::runtime_error("Incorrect magic number: " + filename);
    }
    return result;
}

#ifdef MYRICUBE_WINDOWS
// todo
#else
void* map_file_impl(const filename_string& filename, size_t sz, int* flags)
{
    const char* doing = "";
    int code = 0;
    auto check_code = [&]
    {
        if (code != 0) {
            std::string msg = doing;
            msg += ' ';
            msg += filename;
            msg += ": ";
            msg += strerror(errno);
            throw std::runtime_error(msg);
        }
    };

    // Try to open the file if it exists.
    // int open_flags = *flags & readonly_flag ? O_RDONLY : O_RDWR;
    int open_flags = O_RDWR;
    int fd = open(filename.c_str(), open_flags);

    // Cheap RAII for fd.
    struct Closer {
        int& fd_ref;
        ~Closer() { if (fd_ref >= 0) close(fd_ref); }
    };
    Closer closer = { fd };

    // Depending on create_flag, either return null or create the file
    // if it doesn't exist.
    if (fd < 0) {
        doing = "Opening";
        if (errno != ENOENT) check_code();
        if (!(*flags & create_flag)) return nullptr;

        // Create file now.
        doing = "Creating";
        fd = open(filename.c_str(), O_RDWR | O_CREAT, 0777);
        *flags |= file_created_flag;
        code = fd < 0 ? -1 : 0;
        check_code();

        doing = "Resizing created";
        code = ftruncate(fd, sz);
        check_code();
    }

    // Check the file size is expected, or resize from 0 if needed
    // (and allowed by the create_flag).
    struct stat stats;
    doing = "Checking file size of";
    code = fstat(fd, &stats);
    check_code();

    if (stats.st_size == 0) {
        if (*flags & create_flag) {
            doing = "Resizing";
            *flags |= file_created_flag;
            code = ftruncate(fd, sz);
            check_code();
        }
        else {
            return nullptr;
        }
    }
    else if (stats.st_size != (off_t)sz) {
        throw std::runtime_error(filename +
            " incorrect size (consider removing file manually)");
    }

    // Finally can access the mapping.
    auto prot = PROT_READ | PROT_WRITE;
    // if (!(*flags & readonly_flag)) prot |= PROT_WRITE;
    void* mapping = mmap(nullptr, sz, prot, MAP_SHARED, fd, 0);

    if (mapping == nullptr) {
        doing = "Memory mapping";
        code = -1;
        check_code();
    }
    return mapping;
}
#endif

} // end namespace myricube
