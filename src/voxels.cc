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

#ifdef MYRICUBE_WINDOWS
// todo
#else
#include <arpa/inet.h> /* htonl needed for endian hack */
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

namespace myricube {

// In memory file system entry.
struct InMemoryFile
{
    void* ptr;
    size_t size;
};

// Map from filenames to in-memory files. Nothing is removed from this.
static std::unordered_map<filename_string, InMemoryFile> in_memory_filesystem;

// Set of pointers corresponding to the (start of) in-memory
// files. This is needed because real files need to be unmapped
// (multiple views of the same file use different virtual addresses,
// so this still works with multiple views) while in-memory files
// shall not (this would discard the contents of the in-memory file,
// which is not what I want!)
//
// TODO: Use lower-order pointer bit tricks to remove this???
static std::unordered_set<uintptr_t> in_memory_file_ptrs;

// Mutex protecting above two.
static std::mutex in_memory_filesystem_mutex;



// Given the name of a file (possibly an in-memory temp file), mmap it
// and return the pointer. If the file doesn't exist or has 0 size,
// either return nullptr or create it depending on flags (defined
// below). If created, the T constructor is run in-place on the
// file. THIS CONSTRUCTOR IS EXPECTED TO NOT THROW AN EXCEPTION.
//
// Throw an exception on failure to open for any reason other than
// file not existing. An exception is also thrown if the file size is
// nonzero and doesn't match that of type T (for now I don't
// anticipate mapping arrays).
template <typename T> T* map_file(const filename_string& filename, int* flags);
template <typename T> T* map_file(const filename_string& filename, int flags)
{
#ifndef MYRICUBE_WINDOWS
    // fprintf(stderr, "Mapping %s\n", filename.c_str());
#endif
    auto ptr = map_file<T>(filename, &flags);
    if (!ptr) fprintf(stderr, "map_file returned nullptr.\n");
    return ptr;
}

// Forget about readonly_flag for now because of mutable atomic counter.
// constexpr int readonly_flag = 1; // Added automatically if const ptr requested.

constexpr int create_flag = 2;   // If set, create file if needed.
constexpr int file_created_flag = 4; // Set by the map_file function iff
                                     // the file needed to be created.

void* map_disk_file_impl(const filename_string&, size_t sz, int* flags);
void* map_mem_file_impl(const filename_string&, size_t sz, int* flags);

template <typename T> T* map_file(const filename_string& filename, int* flags)
{
    // if (std::is_const_v<T>) *flags |= readonly_flag;
    void* ptr = starts_with_in_memory_prefix(filename)
              ? map_mem_file_impl(filename, sizeof(T), flags)
              : map_disk_file_impl(filename, sizeof(T), flags);
    if (*flags & file_created_flag) {
        new(ptr) T;
    }
    // TODO think about potential race / interruption condition:
    // create file in separate location and move-in when ready?
    return static_cast<T*>(ptr);
}

#ifdef MYRICUBE_WINDOWS
template <typename T> void unmap_if_disk_file(const T* ptr)
{
    // Only unmap actual on-disk files.
    std::lock_guard guard(in_memory_filesystem_mutex);
    if (in_memory_file_ptrs.count(uintptr_t(ptr))) return;

    throw std::runtime_error("Todo: implement unmap_if_disk_file");
}
#else
template <typename T> void unmap_if_disk_file(const T* ptr)
{
    // fprintf(stderr, "unmap_file of size %ld\n", long(sizeof(T)));

    // Only unmap actual on-disk files.
    std::lock_guard guard(in_memory_filesystem_mutex);
    if (in_memory_file_ptrs.count(uintptr_t(ptr))) return;

    auto code = munmap((void*)ptr, sizeof(T));
    assert(code == 0); // Seems like munmap only fails due to logic errors.
}
#endif

void unmap_mut_bin_chunk_group(BinChunkGroup* ptr)
{
    unmap_if_disk_file(ptr);
}

void unmap_bin_chunk_group(const BinChunkGroup* ptr)
{
    unmap_if_disk_file(ptr);
}

void unmap_bin_group_bitfield(BinGroupBitfield* ptr)
{
    unmap_if_disk_file(ptr);
}

WorldHandle::WorldHandle(const filename_string& arg)
{
    auto on_bad_filename = [&arg]
    {
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
        unmap_if_disk_file<BinWorld>);

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

#ifndef MYRICUBE_WINDOWS
    thread_local long endian_fix_enabled;
    if (!endian_fix_enabled) {
        const char* env = getenv("myricube_endian_fix");
        if (env != nullptr) {
            char* endptr;
            endian_fix_enabled = strtol(env, &endptr, 10);
            if (*endptr != 0) {
                throw std::runtime_error("Could not parse myricube_endian_fix "
                    "environment variable as integer.");
            }
        }
    }

    if (!endian_fix_enabled) {
        throw std::runtime_error("Incorrect endianness: " + filename
            + "\nrun with environment variable myricube_endian_fix=1"
              "\nto fix in-place (backup first!)");
    }

    for (int zH = 0; zH < edge_chunks; ++zH) {
    for (int yH = 0; yH < edge_chunks; ++yH) {
    for (int xH = 0; xH < edge_chunks; ++xH) {
        BinChunk& chunk = group.chunk_array[zH][yH][xH];

        for (int zL = 0; zL < chunk_size; ++zL) {
        for (int yL = 0; yL < chunk_size; ++yL) {
        for (int xL = 0; xL < chunk_size; ++xL) {
            uint32_t* p_voxel = &chunk.voxel_array[zL][yL][xL];
            *p_voxel = htonl(*p_voxel);
            // Just get rid of maybe_fix_endian on Windows.
        }
        }
        }
    }
    }
    }

    group.magic_number |= new_endian_magic;
    fprintf(stderr, "Fixed endianness of %s\n", filename.c_str());
    return true;
#else
    fwprintf(stderr,
        L"Incorrect endianness: %ls"
        "\n(cannot fix on Windows, call me and complain if you need this).",
        filename.c_str()
    );
    throw std::runtime_error("Incorrect endianness");
#endif
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
    auto result = UPtrChunkGroup(map_file<const BinChunkGroup>(filename, 0));
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



#ifdef MYRICUBE_WINDOWS

// Map an actual on-disk file.
void* map_disk_file_impl(const filename_string& filename, size_t sz, int* flags)
{
    throw std::runtime_error("TODO: Implement map_disk_file_impl on Windows");
}

// Map an ephemeral in-memory file.
void* map_mem_file_impl(const filename_string& filename, size_t sz, int* flags)
{
    assert(starts_with_in_memory_prefix(filename));

    std::lock_guard lock(in_memory_filesystem_mutex);

    auto iter = in_memory_filesystem.find(filename);

    // Such an in-memory file exists: check size and return.
    if (iter != in_memory_filesystem.end()) {
        InMemoryFile in_memory_file = iter->second;
        if (in_memory_file.size != sz) {
            fprintf(stderr, "%ls: Incorrect size\n", filename.c_str());
            throw std::runtime_error("In memory file: incorrect size");
        }
        assert(in_memory_file.ptr);
        return in_memory_file.ptr;
    }

    // No such in-memory file exists: either create or return nullptr.
    if ((*flags & create_flag) == 0) {
        return nullptr;
    }
    void* mapping = nullptr;
    try {
        mapping = malloc(sz);
        in_memory_file_ptrs.emplace(uintptr_t(mapping));

        if (mapping == nullptr) throw std::bad_alloc();

        InMemoryFile in_memory_file { mapping, sz };
        in_memory_filesystem.emplace(filename, in_memory_file);
        *flags |= file_created_flag;
    }
    catch (...) {
        if (mapping) {
            free(mapping);
            in_memory_file_ptrs.erase(uintptr_t(mapping));
        }
        throw;
    }
    return mapping;
}

#else
// Map an actual on-disk file.
void* map_disk_file_impl(const filename_string& filename, size_t sz, int* flags)
{
    assert(!starts_with_in_memory_prefix(filename));

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
    doing = "Opening";
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
        // ...but first check that the failure was actually that the
        // file didn't exist.
        if (errno != ENOENT) {
            code = fd;
            check_code();
        }
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

// Map an ephemeral in-memory file.
void* map_mem_file_impl(const filename_string& filename, size_t sz, int* flags)
{
    assert(starts_with_in_memory_prefix(filename));

    std::lock_guard lock(in_memory_filesystem_mutex);

    auto iter = in_memory_filesystem.find(filename);

    // Such an in-memory file exists: check size and return.
    if (iter != in_memory_filesystem.end()) {
        InMemoryFile in_memory_file = iter->second;
        if (in_memory_file.size != sz) {
            throw std::runtime_error(filename + " (in memory) incorrect size");
        }
        assert(in_memory_file.ptr);
        return in_memory_file.ptr;
    }

    // No such in-memory file exists: either create or return nullptr.
    if ((*flags & create_flag) == 0) {
        return nullptr;
    }
    void* mapping = nullptr;
    try {
        auto prot = PROT_READ | PROT_WRITE;
        auto mmap_flags = MAP_PRIVATE | MAP_ANONYMOUS;
        mapping = mmap(nullptr, sz, prot, mmap_flags, -1, 0);
        in_memory_file_ptrs.emplace(uintptr_t(mapping));

        if (mapping == nullptr) throw std::bad_alloc();

        InMemoryFile in_memory_file { mapping, sz };
        in_memory_filesystem.emplace(filename, in_memory_file);
        *flags |= file_created_flag;
    }
    catch (...) {
        if (mapping) {
            munmap(mapping, sz);
            in_memory_file_ptrs.erase(uintptr_t(mapping));
        }
        throw;
    }
    return mapping;
}
#endif

} // end namespace myricube
