// Windows and Linux implementations of map.hh

#include "map.hh"

#include <mutex>
#include <unordered_map>
#include <unordered_set>

namespace myricube {

// In memory file system entry.
struct InMemoryFile
{
    void* ptr;
    size_t size;
};

// Map from filenames to in-memory files. Nothing is removed from this.
static std::unordered_map<filename_string, InMemoryFile> in_memory_filesystem;

#ifndef MYRICUBE_WINDOWS
// Map from mappings to their size.
static std::unordered_map<uintptr_t, size_t> mmap_size_map;
#endif

// Set of pointers corresponding to the (start of) in-memory
// files. This is needed because real files need to be unmapped
// (multiple views of the same file use different virtual addresses,
// so this still works with multiple views) while in-memory files
// shall not (this would discard the contents of the in-memory file,
// which is not what I want!)
//
// TODO: Use lower-order pointer bit tricks to remove this???
static std::unordered_set<uintptr_t> in_memory_file_ptrs;

// Mutex protecting above.
static std::mutex map_mutex;

#ifdef MYRICUBE_WINDOWS
void unmap(const void* ptr)
{
    // Only unmap actual on-disk files.
    std::lock_guard guard(map_mutex);
    if (in_memory_file_ptrs.count(uintptr_t(ptr))) return;

    bool okay = UnmapViewOfFile(ptr);
    if (!okay) {
        fprintf(stderr, "UnmapViewOfFile: %i\n", int(GetLastError()));
        throw std::runtime_error("UnmapViewOfFile failed");
    }
}
#else
void unmap(const void* ptr)
{
    // Only unmap actual on-disk files.
    std::lock_guard guard(map_mutex);
    if (in_memory_file_ptrs.count(uintptr_t(ptr))) return;

    auto it = mmap_size_map.find(uintptr_t(ptr));
    assert(it != mmap_size_map.end());
    auto code = munmap((void*)ptr, it->second);
    assert(code == 0); // Seems like munmap only fails due to logic errors.
}
#endif

#ifdef MYRICUBE_WINDOWS

// Map an actual on-disk file.
void* map_disk_file_impl(
    const filename_string& filename,
    size_t sz, int* flags, size_t* p_array_size)
{
    assert(!starts_with_in_memory_prefix(filename));
    bool readonly = *flags & readonly_flag;

    const char* doing = "";
    auto check_last_error = [&]
    {
        DWORD error = GetLastError();
        if (error) {
            wchar_t buf[256];
            FormatMessageW(
                FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                 NULL,
                 error,
                 MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                 buf,
                 (sizeof(buf) / sizeof(wchar_t)),
                 NULL);

            std::wstring msg = filename;
            msg += L": ";
            msg += buf;
            fprintf(stderr, "%ls\n", msg.c_str());
            throw std::runtime_error(
                std::string("map_disk_file_impl failed on ") + doing);
        }
    };

    // Try to open the file if it exists.
    doing = "Opening";
    HANDLE handle = CreateFileW(
        filename.c_str(),
        readonly ? GENERIC_READ : GENERIC_READ | GENERIC_WRITE,
        FILE_SHARE_WRITE | FILE_SHARE_READ | FILE_SHARE_DELETE,
        NULL,
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL,
        NULL);

    // fprintf(stderr, "File=%ls Handle=%p\n", filename.c_str(), handle);

    // Cheap RAII for handle.
    struct Closer {
        HANDLE& handle_ref;
        ~Closer()
        {
            if (handle_ref != INVALID_HANDLE_VALUE) CloseHandle(handle_ref);
        }
    };
    Closer closer = { handle };

    // Depending on create_flag, either return null or create the file
    // if it doesn't exist.
    if (handle == INVALID_HANDLE_VALUE) {
        // ...but first check that the failure was actually that the
        // file didn't exist.
        if (GetLastError() != ERROR_FILE_NOT_FOUND) {
            check_last_error();
        }
        if (!(*flags & create_flag)) return nullptr;

        // Create file now.
        doing = "Creating";
        SetLastError(0);
        handle = CreateFileW(
            filename.c_str(),
            readonly ? GENERIC_READ : GENERIC_READ | GENERIC_WRITE,
            FILE_SHARE_WRITE | FILE_SHARE_READ | FILE_SHARE_DELETE,
            NULL,
            CREATE_NEW,
            FILE_ATTRIBUTE_NORMAL,
            NULL);
        *flags |= file_created_flag;
        check_last_error();

        // Resize the file to the correct size.
        SetLastError(0);
        doing = "Setting file pointer";
        LONG low = LONG(sz);
        LONG high = LONG(sz >> 32);
        SetFilePointer(handle, low, &high, FILE_BEGIN);
        check_last_error();

        doing = "Resizing file";
        SetEndOfFile(handle);
        check_last_error();
    }

    // Check the file size is expected, or resize from 0 if needed
    // (and allowed by the create_flag).
    doing = "Checking file size of";
    size_t file_size;
    {
        DWORD low;
        DWORD high;
        low = GetFileSize(handle, &high);
        check_last_error();
        file_size = size_t(low) | size_t(high) << 32;
    }

    size_t N = file_size / sz;
    if (file_size == 0) {
        if (*flags & create_flag) {
            *flags |= file_created_flag;

            // Resize the file to the correct size.  This might be
            // redundant on Windows (CreateFileMapping resizes
            // implicitly?), not sure.
            SetLastError(0);
            doing = "Setting existing file pointer";
            LONG low = LONG(sz);
            LONG high = LONG(sz >> 32);
            SetFilePointer(handle, low, &high, FILE_BEGIN);
            check_last_error();

            doing = "Resizing existing file";
            SetEndOfFile(handle);
            check_last_error();
        }
        else {
            return nullptr;
        }
    }
    else if (N * sz != file_size or (!p_array_size and N != 1)) {
        fprintf(stderr, "Incorrect file size: %ls\n", filename.c_str());
        throw std::runtime_error(
            "Incorrect file size (consider removing file manually)");
    }
    if (p_array_size) *p_array_size = N;

    // Now I need to create a "file-mapping object". There seems to be
    // no Linux equivalent so let's hope I did this correctly...
    doing = "Creating Windows File Mapping Object";
    HANDLE mapobj = CreateFileMappingA(
        handle, // mabye I should have given it a better name...
        NULL,
        readonly ? PAGE_READ : PAGE_READWRITE,
        DWORD(file_size >> 32), // These are now swapped in order from before!!!
        DWORD(file_size),       // Now low is second instead of first WTF
        NULL);
    Closer mapobj_closer { mapobj }; // More cheap RAII.
    check_last_error();

    // Finally can access the mapping.
    doing = "Memory mapping";
    void* mapping = MapViewOfFile(
        mapobj,
        FILE_MAP_ALL_ACCESS,
        0, 0,
        file_size); // OMG now I can just use size_t directly!!!
    if (mapping == nullptr) check_last_error();
    return mapping;

    // NOTE: At this point the Closer objects destroy the mapping
    // object and the file handle. Some of the docs I read say I
    // shouldn't do this, but this seems based on the assumption that
    // I want exclusive file access, which I don't. According to this
    // stack overflow answer:
    //
    // https://stackoverflow.com/questions/36495158/c-windows-api-close-file-handle-before-unmapviewoffile
    //
    // ...it should be okay to close those handles now.
}

// Map an ephemeral in-memory file.
void* map_mem_file_impl(
    const filename_string& filename,
    size_t sz, int* flags, size_t* p_array_size)
{
    assert(starts_with_in_memory_prefix(filename));

    std::lock_guard lock(map_mutex);

    auto iter = in_memory_filesystem.find(filename);

    // Such an in-memory file exists: check size and return.
    if (iter != in_memory_filesystem.end()) {
        size_t N = file_size / sz;
        InMemoryFile in_memory_file = iter->second;
        if (in_memory_file.size != N * sz or (!p_array_size and N != 1)) {
            fprintf(stderr, "%ls: Incorrect size\n", filename.c_str());
            throw std::runtime_error("In memory file: incorrect size");
        }
        if (p_array_size) *p_array_size = N;
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
    if (p_array_size) *p_array_size = 1;
    return mapping;
}


// Linux implementation of memory-mapped files.
#else
// Map an actual on-disk file.
void* map_disk_file_impl(
    const filename_string& filename,
    size_t sz, int* flags, size_t* p_array_size)
{
    assert(!starts_with_in_memory_prefix(filename));
    bool readonly = *flags & readonly_flag;

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
    doing = "Opening";
    int open_flags = readonly ? O_RDONLY : O_RDWR;
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
        if (!(*flags & create_flag)) return nullptr;

        // Create file now.
        doing = "Creating";
        fd = open(filename.c_str(), open_flags | O_CREAT, 0777);
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

    size_t file_size = size_t(stats.st_size);
    size_t N = file_size / sz;
    if (file_size == 0) {
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
    else if (file_size != N * sz or (!p_array_size and N != 1)) {
        throw std::runtime_error(filename +
            " incorrect size (consider removing file manually)");
    }

    // Finally can access the mapping.
    auto prot = readonly ? PROT_READ : PROT_READ | PROT_WRITE;
    void* mapping = mmap(nullptr, file_size, prot, MAP_SHARED, fd, 0);

    if (mapping == nullptr) {
        doing = "Memory mapping";
        code = -1;
        check_code();
    }
    try {
        std::lock_guard lock(map_mutex);
        mmap_size_map[uintptr_t(mapping)] = file_size;
    }
    catch (...) {
        munmap(mapping, file_size);
        throw;
    }
    if (p_array_size) *p_array_size = N;
    return mapping;
}

// Map an ephemeral in-memory file.
void* map_mem_file_impl(
    const filename_string& filename,
    size_t sz, int* flags, size_t* p_array_size)
{
    assert(starts_with_in_memory_prefix(filename));

    std::lock_guard lock(map_mutex);

    auto iter = in_memory_filesystem.find(filename);

    // Such an in-memory file exists: check size and return.
    if (iter != in_memory_filesystem.end()) {
        InMemoryFile in_memory_file = iter->second;
        size_t N = in_memory_file.size / sz;
        if (in_memory_file.size != N * sz or (!p_array_size and N != 1)) {
            throw std::runtime_error(filename + " (in memory) incorrect size");
        }
        if (p_array_size) *p_array_size = N;
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
    if (p_array_size) *p_array_size = 1;
    return mapping;
}
#endif

} // end namespace myricube
