// Cross-platform file mapping utility. "Files" whose names start with
// the in-memory prefix (defined later) are quietly made to refer to
// chunks of allocated memory (lost upon program shutdown) rather than
// actual files.

#ifndef MYRICUBE_MAP_HH_
#define MYRICUBE_MAP_HH_

#include "myricube.hh"

#include <stdio.h>
#include <string>

#ifdef MYRICUBE_WINDOWS
#include <windows.h>

#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#endif

namespace myricube {

// /MEM:/ is the in_memory_prefix. Colons are not allowed in file
// names on Windows and almost certainly not in the name of a
// top-level directory on Unix.
inline const filename_string get_in_memory_prefix()
{
#ifdef MYRICUBE_WINDOWS
    return L"/MEM:/";
#else
    return "/MEM:/";
#endif
}

template <typename String>
inline bool starts_with_in_memory_prefix(const String& arg)
{
    filename_string in_memory_prefix = get_in_memory_prefix();
    if (arg.size() < in_memory_prefix.size()) return false;
    size_t index = 0;
    for (auto c : in_memory_prefix) {
        if (arg[index++] != c) return false;
    }
    return true;
}

inline filename_string add_in_memory_prefix(const filename_string& arg)
{
    if (starts_with_in_memory_prefix(arg)) return arg;

    filename_string result = get_in_memory_prefix();
    result.append(arg);
    return result;
}

#ifdef MYRICUBE_WINDOWS
inline filename_string add_in_memory_prefix(const std::string& arg)
{
    if (starts_with_in_memory_prefix(arg)) {
        filename_string result;
        result.reserve(arg.size());
        for (char c : arg) result.push_back(wchar_t(c));
        return result;
    }

    return filename_concat_c_str(get_in_memory_prefix(), arg.c_str());
}
#endif

// Given the name of a file (possibly an in-memory temp file), mmap it
// and return the pointer. If the file doesn't exist or has 0 size,
// either return nullptr or create it depending on flags (defined
// below). If created, the file is of size equal to T and the T
// constructor is run in-place on the file. THIS CONSTRUCTOR IS
// EXPECTED TO NOT THROW AN EXCEPTION.
//
// If p_array_size is provided, *p_array_size is set to
//    N = (size of file) / (size of T).
//
// Throw an exception on failure to open for any reason other than
// file not existing, or if N is not an integer. If p_array_size is not
// provided, an exception is throw if N is not either 0 or 1.

template <typename T>
T* map_file(
    const filename_string& filename, int* flags, size_t* p_array_size=nullptr);

template <typename T>
T* map_file(
    const filename_string& filename, int flags, size_t* p_array_size=nullptr)
{
    auto ptr = map_file<T>(filename, &flags, p_array_size);
    return ptr;
}

// Flags for map_file
constexpr int readonly_flag = 1;
constexpr int create_flag = 2;   // If set, create file if needed.
constexpr int file_created_flag = 4; // Set by the map_file function iff
                                     // the file needed to be created.

void* map_disk_file_impl(
    const filename_string&, size_t sz, int* flags, size_t* p_array_size);
void* map_mem_file_impl(
    const filename_string&, size_t sz, int* flags, size_t* p_array_size);

template <typename T>
T* map_file(const filename_string& filename, int* flags, size_t* p_array_size)
{
    void* ptr = starts_with_in_memory_prefix(filename)
              ? map_mem_file_impl(filename, sizeof(T), flags, p_array_size)
              : map_disk_file_impl(filename, sizeof(T), flags, p_array_size);
    if (*flags & file_created_flag) {
        assert(ptr != nullptr);
        new(ptr) T;
    }
    // TODO think about potential race / interruption condition:
    // create file in separate location and move-in when ready?
    return static_cast<T*>(ptr);
}

// Unmap previous result of map_file. In-memory files persist (until
// shutdown) even if unmapped; they can be remapped again.
void unmap(const void* ptr);

} // end namespace myricube

#endif /* !MYRICUBE_MAP_HH_ */