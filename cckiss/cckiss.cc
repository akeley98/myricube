// By David Zhao Akeley (dza724`at`gmail.com / akeley98`at`g.ucla.edu) 2019
// Released under dual license: your choice of
//
// 1. Creative Commons CC0 license ("No Rights Reserved")
//    https://creativecommons.org/publicdomain/zero/1.0/legalcode
//
// 2. Do What the Fuck You Want to Public License http://www.wtfpl.net/about/
//
// Summary:
//
// 1. This program takes as its first argument the name of a target
// file to compile. The target file name should be of the form
// "cckiss/[source file].[s|o]", where [source file] is the path to
// the C, C++, or GLSL file that compiles to the target. Anything with
// a .glsl extension is assumed GLSL, .c extension is assumed to be C;
// all others C++. Further arguments control:
//
//     a. The C/C++ compiler name (prefix with ".cckiss.CXX", or first
//     argument after target filename).
//
//     b. Arguments for invoking the compiler as a preprocessor
//     (prefix with ".cckiss.CPPFLAGS")
//
//     c. Arguments for invoking the C/C++ compiler (prefix with
//     ".cckiss.CXXFLAGS", even if we are using a C compiler). Note:
//     '-c' or '-S' need not be included: it will be added
//     automatically depending on whether the target file is an object
//     (.o) or assembly (.s) file.
//
//     d. The GLSL compiler (prefix with .cckiss.GLSLC), if needed.
//
//     e. Arguments for the GLSL compiler (prefix with .cckiss.GLSLARGS).
//
// 2. If the target need not be recompiled, the program terminates. Otherwise,
//
// 3. The source file is preprocessed, with the preprocessor output stored in
// "cckiss/[source file].[i|ii]" (extension differs for C/C++).
//
// 4. The preprocessed file is scanned for preprocessor directives of
// the form "#{line} [line number] "[included file name]"". This gives
// us a list of all dependency files, which is stored in
// "cckiss/[source file]-deps.txt".  Exactly one newline separates
// each file listed.
//
// 5. Future invokations of cckiss can check that list of dependency files,
// and skip recompilation if all dependency files' modification times are
// before the timestamp of the target file. (However, recompilation is done
// unconditionally if 'B' is in the environment variable 'MAKEFLAGS'. This
// is for compatibility with make's -B argument).

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <fstream>
#include <set>
#include <stdio.h>
#include <stdexcept>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <string>
#include <string.h>
#include <unistd.h>
#include <utility>
#include <vector>

#define RED "\x1b[1m\x1b[31m"
#define BLUE "\x1b[1m\x1b[34m"
#define PURPLE "\x1b[35m\x1b[1m"
#define CYAN "\x1b[36m\x1b[1m"
#define END "\x1b[0m"

namespace {

std::string cxxflags_delim = ".cckiss.CXXFLAGS";
std::string cppflags_delim = ".cckiss.CPPFLAGS";
std::string cxx_delim = ".cckiss.CXX";
std::string glslc_delim = ".cckiss.GLSLC";
std::string glslargs_delim = ".cckiss.GLSLARGS";

enum class FileType { zero, c, cxx, glsl };

struct Args
{
    // Name of compiled file.
    std::string target_file_name;

    // Name of the source file compiled to target file.
    std::string source_file_name;

    // Name of the file preprocessor output for the source file is written to.
    std::string preprocessed_file_name;

    // Name of the dependencies file scanned from the preprocessed file.
    std::string deps_file_name;

    // Makefile CC or CXX
    std::string cxx;

    // Name of a dependency file whose change triggered recompilation.
    // Filled in later; if empty and (re)compiling is done, the compilation
    // was caused by some other reason (e.g. target never was compiled)
    std::string changed_dependency_name;

    // Makefile CPPFLAGS
    std::vector<std::string> cppflags;

    // Makefile CFLAGS or CXXFLAGS
    std::vector<std::string> cxxflags;

    // Makefile GLSLC
    std::string glslc;

    // Makefile GLSLARGS
    std::vector<std::string> glslargs;

    FileType source_file_type = FileType::zero;

    bool verbose = false;
    bool always_make = false;
};

// Arguably evil, but allow timespecs to be compared like ordinary
// darn numbers.
inline bool operator< (struct timespec left, struct timespec right) {
    return left.tv_sec < right.tv_sec ||
        (left.tv_sec == right.tv_sec && left.tv_nsec < right.tv_nsec);
}
inline bool operator> (struct timespec left, struct timespec right) {
    return right < left;
}
inline bool operator<= (struct timespec left, struct timespec right) {
    return !(left > right);
}
inline bool operator>= (struct timespec left, struct timespec right) {
    return !(left < right);
}

// If str starts with the string pointed to by prefix, return a pointer
// to the remaining part of str. Otherwise return null.
const char* skip_prefix(const char* str, const char* prefix) {
    while (*prefix != 0) {
        if (*str++ != *prefix++) return 0;
    }
    return str;
}

// For the named file, return
//
// 1. Whether the file exists, and
// 2. If so, the modification time, through *timespec.
//
// Terminates if the file exists but could not be accessed.
bool file_exists_mtim(const std::string& name, struct timespec *timespec)
{
    struct stat statbuf;
    int return_code = 1;
    constexpr int max_tries = 10000;
    for (int tries = 0; return_code != 0; ++tries) {
        return_code = stat(name.c_str(), &statbuf);
        if (return_code == 0) break;

        // File doesn't exist:
        if (errno == ENOENT) return false;

        // Unix stupidity:
        if (errno == EINTR && tries < max_tries) continue;

        const char* msg = strerror(errno);
        printf("Could not access \"%s\": %s\n", name.c_str(), msg);
        fprintf(stderr, "Could not access \"%s\": %s\n", name.c_str(), msg);
        exit(1);
    }
    *timespec = statbuf.st_mtim;
    return true;
}

// Given the name of a target file (compiled object/assembly file),
// return the name of the corresponding source file. Example:
// "cckiss/pear/util.cc.s" -> "pear/util.cc".
std::string source_file_name_from_target(const std::string& compiled_file_name)
{
    // Need to recover the name of the source file for the given
    // target (compiled) file. Expect compiled_file_name =
    // "cckiss/{source/file}.[s|o]"
    std::string source_file_name;
    size_t sz;
    const char* skipped_prefix = skip_prefix(
        compiled_file_name.c_str(),
        "cckiss/");
    if (!skipped_prefix) goto bad_file_name;

    // Now source_file_name = {source/file}.[s|o] Remove the .s or .o
    source_file_name = std::string(skipped_prefix);
    sz = source_file_name.size();
    if (sz <= 2) goto bad_file_name;
    if (source_file_name[sz-1] != 's' && source_file_name[sz-1] != 'o') {
        goto bad_file_name;
    }
    if (source_file_name[sz-2] != '.') goto bad_file_name;
    source_file_name.pop_back();
    source_file_name.pop_back();
    return source_file_name;

bad_file_name:
    printf(
        "Expected target (compiled) file name to be of the form "
        "cckiss/{path/to/source/file}.[s|o]. Got: \"%s\".\n",
        compiled_file_name.c_str());
    fprintf(
        stderr,
        "Expected target (compiled) file name to be of the form "
        "cckiss/{path/to/source/file}.[s|o]. Got: \"%s\".\n",
        compiled_file_name.c_str());
    exit(1);
}

// Given the name of a source file, return the expected filename of
// its dependencies file (which was created by cckiss based on scanning
// the source file post-preprocessing)
std::string deps_file_name_for_source(
    const std::string& source_file_name)
{
    return "cckiss/" + source_file_name + "-deps.txt";
}

// Given the name and type of a source file, return the expected
// filename of its corresponding preprocessed file.
std::string preprocessed_file_name_for_source(
    const std::string& source_file_name, FileType file_type)
{
    // C files must become .i files -- the only file extension I know
    // of for C is lowercase '.c'.

    // GLSL files remain .glsl; only difference is the preprocessed file
    // is in the cckiss directory.

    // All others are assumed to be C++, with .ii extension once
    // preprocessed.
    const char* extension = "";
    switch (file_type) {
        default: assert(0); break;
        case FileType::c: extension = ".i"; break;
        case FileType::cxx: extension = ".ii"; break;
        case FileType::glsl: extension = ""; break; // has .glsl already.
    }
    return "cckiss/" + source_file_name + extension;
}

// The cckiss/ directory contains a replica of the source directory
// tree (e.g. if we are compiling foo/bar/barf.c, the cckiss files
// generated for it -- barf.c.i, barf.c-deps.txt, barf.c.s -- are
// stored in cckiss/foo/bar/).
//
// This function ensures that the needed directories in cckiss/ are
// created for the named source file.
void make_directory_tree_for_source(const std::string& source_file_name)
{
    // Slowly copy `source_file_name` to `directory_to_create`,
    // pausing at each '/' to create the directory.
    std::string directory_to_create = "cckiss/";
    const char* pathname = "";
    for (char c : source_file_name) {
        directory_to_create += c;
        if (c == '/') {
            pathname = directory_to_create.c_str();
            int return_code = mkdir(pathname, 0777);
            bool eexist = (errno == EEXIST);
            if (return_code < 0 && !eexist) {
                goto fail;
            }
            if (eexist) {
            retry:
                struct stat statbuf;
                return_code = stat(pathname, &statbuf);
                if (return_code < 0) {
                    if (errno == EINTR) goto retry;
                    else goto fail;
                }
                if (!S_ISDIR(statbuf.st_mode)) {
                    errno = EISDIR;
                    goto fail;
                }
            } else {
                printf("# Created directory \"%s\".\n", pathname);
            }
        }
    }
    return;
fail:
    const char* msg = strerror(errno);
    printf("Failed to create \"%s\": %s\n", pathname, msg);
    fprintf(stderr, "Failed to create \"%s\": %s\n", pathname, msg);
    exit(1);
}

// Read the dependency file correspoding to the arg's source file, and
// write to the vector the list of dependency file names within
// (newline separated in deps file). Iff the dependency file doesn't
// yet exist, return false. If the dependency file cannot be read for
// any other reason, terminate.
bool try_parse_deps_file(Args& args, std::vector<std::string>* out)
{
    auto& deps = *out;
    deps.clear();
    auto& source_file_name = args.source_file_name;
    auto& deps_file_name = args.deps_file_name;

    std::ifstream stream(deps_file_name);
    if (!stream.is_open()) {
        struct stat statbuf;
        int return_code = stat(deps_file_name.c_str(), &statbuf);

        // Check "impossible" error: fstream can't open the file but
        // stat can (race condition possible I suppose).
        if (return_code == 0) {
            printf(
                "Could not read \"%s\" for unknown reason.\n",
                source_file_name.c_str());
            fprintf(
                stderr,
                "Could not read \"%s\" for unknown reason.\n",
                source_file_name.c_str());
            exit(1);
        }

        // Failed to open deps file because it doesn't exist: okay,
        // just report it to the caller.
        if (errno == ENOENT) {
            return false;
        }

        // Otherwise, something else isn't right and needs to be reported.
        const char* reason = strerror(errno);
        if (return_code == 0) {
            printf(
                "Could not read \"%s\": %s\n",
                source_file_name.c_str(),
                reason);
            fprintf(
                stderr,
                "Could not read \"%s\": %s\n",
                source_file_name.c_str(),
                reason);
            exit(1);
        }
    }

    std::string dep_file_name;

    while (1) {
        std::getline(stream, dep_file_name, '\n');
        if (stream.eof()) return true;
        if (stream.bad()) {
            fprintf(
                stderr,
                "Internal error: ifstream badbit set for \"%s\".\n",
                deps_file_name.c_str());
            exit(1);
        }
        if (stream.fail()) {
            fprintf(
                stderr,
                "Internal error: ifstream failbit set for \"%s\".\n",
                deps_file_name.c_str());
            exit(1);
        }
        deps.push_back(std::move(dep_file_name));
    }
}

// Look at args for the name of the target file (compiled
// object/assembly file in cckiss/), and return true if we should
// recompile it. This is the case if any of the following hold:
//
// 1. Target file doesn't yet exist.
//
// 2. The target's corresponding source file was modified.
//
// 3. One of the dependencies listed in the source file's
//    corresponding dependency file was modified.
//
// 4. Said dependency file does not exist.
bool should_recompile_target_file(Args& args)
{
    auto& compiled_file_name = args.target_file_name;
    auto& source_file_name = args.source_file_name;

    struct timespec target_mtim;
    bool compiled_file_exists = false;
    struct timespec dep_mtim;

    // Check that the required source file exists, and get its
    // modification time. We do this first to provide a warning if the
    // source file doesn't exist (instead of failing mysteriously in
    // the preprocessing stage).
    bool source_exists = file_exists_mtim(source_file_name, &dep_mtim);
    if (!source_exists) {
        printf(
            "# Missing needed source file \"%s\".\n",
            source_file_name.c_str());
        fprintf(
            stderr,
            "Missing needed source file " RED "%s\n" END,
            source_file_name.c_str());
        bool has_file_extension = false;
        for (char c : source_file_name) has_file_extension |= (c == '.');

        if (!has_file_extension) {
            fprintf(stderr,
                "Did you mean to build cckiss/%s" RED ".c" END ".%c or "
                "something instead of %s?\n",
                source_file_name.c_str(),
                compiled_file_name.back(),
                compiled_file_name.c_str());
        }
        exit(1);
    }

    // Extract the modification time of the target file (if it
    // exists). Later, compare this modification time with the
    // modification times of the dependency files.
    compiled_file_exists =
        file_exists_mtim(compiled_file_name, &target_mtim);

    // If the target file doesn't exist yet, we obviously must compile.
    if (!compiled_file_exists) {
        printf("# New target file \"%s\".\n", compiled_file_name.c_str());
        return true;
    }
    else {
        // Otherwise, compare the dependency modification times with the
        // target file modification time.

        // Check source file (dep_mtim should have been filled in earlier).
        if (dep_mtim >= target_mtim) {
            printf(
                "# \"%s\" modified, needed by \"%s\".\n",
                source_file_name.c_str(),
                compiled_file_name.c_str());
            args.changed_dependency_name = source_file_name;
            return true;
        }

        // Check dependency files. If we don't have the list of
        // dependencies, we have to recompile to create it.
        std::vector<std::string> deps_file_names;
        bool has_deps_file = try_parse_deps_file(args, &deps_file_names);
        if (!has_deps_file) return true;

        for (const std::string& dep_file_name : deps_file_names) {
            // Skip blank lines, lines starting with '<', and ending with '/'.
            // This is random cruft that compilers stick in. TODO
            if (dep_file_name.size() == 0) continue;
            if (dep_file_name.front() == '<' || dep_file_name.back() == '/') {
                if (args.verbose) {
                    printf("# Ignoring file name \"%s\"\n",
                        dep_file_name.c_str());
                }
                continue;
            }

            bool dep_exists = file_exists_mtim(dep_file_name, &dep_mtim);
            if (!dep_exists) {
                printf("# Ignoring missing \"%s\", "
                    "possibly needed by \"%s\".\n",
                    dep_file_name.c_str(),
                    compiled_file_name.c_str());
            }

            if (dep_mtim >= target_mtim) {
                printf("# Dependency \"%s\" modified; needed by \"%s\".\n",
                    dep_file_name.c_str(),
                    compiled_file_name.c_str());
                args.changed_dependency_name = dep_file_name;
                return true;
            }
        }

        // No dependency changes detected at this point.
        return false;
    }
}

int preprocess_source_to_fd(Args& args);
void make_deps_file_from_fd(Args& args, int fd);

// Look at args for the name of the source file, preprocess it
// (storing the preprocessed file under the corresponding directory in
// "cckiss/"), and scan the preprocessed file for dependency file
// names (storing the deps file in "cckiss/").
void preprocess_and_make_deps_file(Args& args)
{
    make_directory_tree_for_source(args.source_file_name);
    int fd = preprocess_source_to_fd(args);
    try {
        make_deps_file_from_fd(args, fd);
    }
    catch (...) {
        close(fd);
        throw;
    }

    close(fd);
}

// Given the command line args, name of a source file, and the name of
// the preprocessed file, preprocess the named source file. Returns
// the file descriptor of the preprocessed file. Assumes that CXX or
// GLSLC can be run as a preprocessor with '-E' argument, and that it
// writes the output to stdout (fd 1).
int preprocess_source_to_fd(Args& args)
{
    // First, open the preprocessor output file.
    const char* pathname = args.preprocessed_file_name.c_str();
retry:
    int fd = open(pathname, O_RDWR | O_CREAT | O_TRUNC, 0666);
    if (fd < 0) {
        if (errno == EINTR) goto retry;
        const char* msg = strerror(errno);
        printf("Could not open \"%s\": %s.\n", pathname, msg);
        fprintf(stderr, "Could not open \"%s\": %s.\n", pathname, msg);
        exit(1);
    }

    // Use C/C++ compiler for C/C++ files, glsl compiler for glsl.
    assert(args.source_file_type != FileType::zero);
    const char* compiler =
        args.source_file_type == FileType::glsl ?
        args.glslc.c_str() : args.cxx.c_str();
    const std::vector<std::string>& cxxflags =
        args.source_file_type == FileType::glsl ?
        std::vector<std::string>{} : args.cxxflags;
    const std::vector<std::string>& cppflags =
        args.source_file_type == FileType::glsl ?
        args.glslargs : args.cppflags;

    // Convert the std::string args into a list of char pointers.
    std::vector<const char*> argv;
    argv.push_back(compiler);
    for (const auto& arg_string : cxxflags) {
        // Vitally important that arg_string is by-reference here.
        argv.push_back(arg_string.c_str());
    }
    for (const auto& arg_string : cppflags) {
        argv.push_back(arg_string.c_str());
    }
    static const char _e[] = "-E";
    argv.push_back(_e);
    argv.push_back(args.source_file_name.c_str());

    const char* fmt = "%s";
    for (const char* arg : argv) {
        printf(fmt, arg);
        fmt = " %s";
    }
    printf("\n");
    argv.push_back(nullptr);

    fflush(stdout);
    fflush(stderr);
    auto pid = fork();

    // Error
    if (pid < 0) {
        const char* msg = strerror(errno);
        printf("fork failed: %s.\n", msg);
        fprintf(stderr, "fork failed: %s.\n", msg);
        exit(1);
    }
    // Child
    else if (pid == 0) {
    retry_dup2:
        auto dup2_code = dup2(fd, 1);
        if (dup2_code < 0) {
            if (errno == EINTR) goto retry_dup2;
            const char* msg = strerror(errno);
            printf("dup2 failed: %s.\n", msg);
            fprintf(stderr, "dup2 failed: %s.\n", msg);
            _exit(1);
        }

        char** argv_ptr = const_cast<char**>(argv.data());
        execvp(compiler, argv_ptr);
        const char* msg = strerror(errno);
        printf("Execute %s failed: %s.\n", compiler, msg);
        fprintf(stderr, "Execute %s failed: %s.\n", compiler, msg);
        _exit(2);
    }
    // Parent
    int wstatus;
    auto waitpid_code = waitpid(pid, &wstatus, 0);
    if (waitpid_code < 0) {
        const char* msg = strerror(errno);
        printf("waitpid failed: %s.\n", msg);
        fprintf(stderr, "waitpid failed: %s.\n", msg);
        exit(1);
    }
    if (wstatus != 0) {
        printf("Preprocessing failed: status %i.\n", wstatus);
        fprintf(stderr, "Preprocessing failed: status %i.\n", wstatus);
        exit(WEXITSTATUS(wstatus) || 1);
    }
    return fd;
}

// Return true iff the given line_text is a string of the form:
//
//     [whitespace]#[whitespace][number] "[filename(any str)]"[crud]
//
// or
//
//     [whitespace]#line[whitespace][number] "[filename(any str)]"[crud]
//
// I could use a regex, but I didn't for now (I don't trust myself to
// use them correctly). If it is of the correct form, write to
// *dependency_file_name the contents of the [filename (any str)]
// space above.
bool interpret_as_file_directive(
    Args& args,
    const std::string& line_text,
    long line_number,
    std::string* dependency_file_name)
{
    const char* ptr = line_text.c_str();
    char c = '\0';

    // Expect [whitespace]#[whitespace] or
    // [whitespace]#line[whitespace]. One or both whitespaces may be
    // empty. The magic do loop skips ptr to point to the first
    // non-whitespace char.
    do { c = *ptr; } while (isspace(c) && ++ptr);
    if ( *ptr++ != '#') return false;
    if (const char* skipped = skip_prefix(ptr, "line")) ptr = skipped;
    do { c = *ptr; } while (isspace(c) && ++ptr);

    // Expect number + optional whitespace, and at least one digit.
    c = *ptr;
    if (!isdigit(c)) return false;
    do { c = *ptr; } while (isdigit(c) && ++ptr);
    do { c = *ptr; } while (isspace(c) && ++ptr);

    // Expect double quote.
    if (*ptr++ != '"') return false;

    // Copy out the string inside the double quotes.
    dependency_file_name->clear();
    while (1) {
        c = *ptr++;
        if (c == '"') break;
        if (c == '\0') goto missing_trailing_quote;
        *dependency_file_name += c;
    }

    // Warn for trailing cruft.
    do { c = *ptr; } while ((isspace(c) || isdigit(c)) && ++ptr);
    if (c != '\0') {
        fprintf(
            stderr,
            PURPLE "Suspicious cruft at %s:%li '%s'\n" END,
            args.preprocessed_file_name.c_str(),
            line_number,
            line_text.c_str());
    }
    return true;

missing_trailing_quote:
    fprintf(
        stderr,
        PURPLE "Ignored %s:%li because of missing '\"' in '%s'\n" END,
        args.preprocessed_file_name.c_str(),
        line_number,
        line_text.c_str());
    return false;
}

// Given the file descriptor of the preprocessed file, scan it for
// dependency files (using preprocessor directives of the form
// '# [number] [filename]'). Write the newline-separated list of
// dependency files to the deps file (named in args).
void make_deps_file_from_fd(Args& args, int fd)
{
    // Rewind the fd.
    auto lseek_code = lseek(fd, 0, SEEK_SET);
    if (lseek_code < 0) {
        const char* msg = strerror(errno);
        printf("lseek failed: %s.\n", msg);
        fprintf(stderr, "lseek failed: %s.\n", msg);
        exit(1);
    }

    // I'm using old-fashioned manually buffered IO due to stubborness
    // (I got my brownie points for using ifstream once, now leave me be).
    char buffer[32768];
    char* next_char = &buffer[0];
    char* endptr = next_char;
    bool eof = false;
    auto get_char = [fd, &buffer, &next_char, &endptr, &eof] () -> char
    {
        if (next_char >= endptr) {
        retry:
            auto bytes_read = read(fd, buffer, sizeof buffer);
            if (bytes_read == 0) {
                eof = true;
                return '\n';
            }
            else if (bytes_read < 0) {
                if (errno == EINTR) goto retry;
                const char* msg = strerror(errno);
                printf("read preprocessed file failed: %s.\n", msg);
                fprintf(stderr, "read preprocessed file failed: %s.\n", msg);
                exit(1);
            }
            else {
                next_char = &buffer[0];
                endptr = &buffer[bytes_read];
            }
        }
        char c = *next_char++;
        return c;
    };

    // Process input line-by-line.
    long line_number = 0;
    std::string line_text;
    auto get_line = [&line_text, &line_number, &get_char, &eof] ()
    {
        line_text.clear();
        ++line_number;
        while (!eof) {
            char c = get_char();
            if (eof || c == '\n') return;
            else line_text += c;
        }
    };

    std::set<std::string> deps_set;
    std::string dependency_file_name;

    do {
        get_line();
        bool was_dependency_line =
            interpret_as_file_directive(
                args, line_text, line_number, &dependency_file_name);
        if (was_dependency_line) {
            auto pair = deps_set.insert(dependency_file_name);
            if (pair.second && args.verbose) {
                printf("    %s:%li\n    \"%s\"\n",
                    args.preprocessed_file_name.c_str(),
                    line_number,
                    dependency_file_name.c_str());
            }
        }
    } while (!eof);

    std::string tmp_deps_file_name = args.deps_file_name + "~";

    // Write all the dependency files to a temp file, then move the
    // temp file to the dependencies file location.
retry_open:
    const char* pathname = tmp_deps_file_name.c_str();
    fd = open(pathname, O_CREAT|O_TRUNC|O_WRONLY, 0666);
    if (fd < 0) {
        if (errno == EINTR) goto retry_open;
        const char* msg = strerror(errno);
        printf("Could not open \"%s\": %s\n", pathname, msg);
        fprintf(stderr, "Could not open \"%s\": %s\n", pathname, msg);
        exit(1);
    }

    for (const auto& dep : deps_set) {
        std::string str_to_write = dep + '\n';
        const char* ptr = str_to_write.c_str();
        ssize_t bytes_left = ssize_t(str_to_write.size());
        ssize_t bytes_written;
    retry_write:
        bytes_written = write(fd, ptr, bytes_left);
        if (bytes_written < 0) {
            const char* msg = strerror(errno);
            printf("Could not write to \"%s\": %s\n", pathname, msg);
            fprintf(stderr, "Could not write to \"%s\": %s\n", pathname, msg);
            exit(1);
        }
        ptr += bytes_written;
        bytes_left -= bytes_written;
        if (bytes_left != 0) goto retry_write;
    }

    auto close_code = close(fd);
    if (close_code < 0) {
        const char* msg = strerror(errno);
        printf("Could not close \"%s\": %s\n", pathname, msg);
        fprintf(stderr, "Could not close \"%s\": %s\n", pathname, msg);
        exit(1);
    }

retry_rename:
    auto rename_code = rename(pathname, args.deps_file_name.c_str());
    if (rename_code < 0) {
        if (errno == EINTR) goto retry_rename;
        const char* msg = strerror(errno);
        printf("Could not rename: %s\n", msg);
        fprintf(stderr, "Could not rename: %s\n", msg);
        exit(1);
    }
}

// Run the glsl compiler, converting the preprocessed glsl file to
// a C array of SPIR-V byte code.
//
// The C file's name is just the glsl preprocessed file name + ".c"
// This also generates an included C file (.inc extension).
//
// The array's name is the original glsl source file, with non
// alphanumeric characters converted to underscores.
//
// We also append a size_t to the source file that is the size of
// the array in bytes. This is "sizeof_" + array variable name.
void compile_glsl_to_c(Args& args)
{
    assert(args.source_file_type == FileType::glsl);

    const char* compiler = args.glslc.c_str();
    std::vector<const char*> argv = { compiler };

    for (const auto& arg : args.glslargs) {
        argv.push_back(arg.c_str()); // Note arg by-reference
    }

    // -V
    const char V[] = "-V";
    argv.push_back(V);

    // --vn [array variable name]
    const char vn[] = "--vn";
    argv.push_back(vn);
    std::string variable_name = args.source_file_name;
    for (char& c : variable_name) if (!isalnum(c)) c = '_';
    argv.push_back(variable_name.c_str());

    // -o [output .inc file name]
    const char o[] = "-o";
    argv.push_back(o);
    std::string inc_file_name = args.preprocessed_file_name + ".inc";
    argv.push_back(inc_file_name.c_str());

    // Input file name
    argv.push_back(args.preprocessed_file_name.c_str());

    // Print command and add terminating null pointer.
    const char* fmt = "%s";
    for (const char* arg : argv) {
        printf(fmt, arg);
        fmt = " %s";
    }
    printf("\n");
    argv.push_back(nullptr);

    auto pid = fork();

    // Error
    if (pid < 0) {
        const char* msg = strerror(errno);
        printf("fork failed: %s.\n", msg);
        fprintf(stderr, "fork failed: %s.\n", msg);
        exit(1);
    }
    // Child
    else if (pid == 0) {
        char** argv_ptr = const_cast<char**>(argv.data());
        execvp(compiler, argv_ptr);
        const char* msg = strerror(errno);
        printf("Execute %s failed: %s.\n", compiler, msg);
        fprintf(stderr, "Execute %s failed: %s.\n", compiler, msg);
        _exit(2);
    }

    // While the compiler runs, output the .c file. This just includes
    // headers for uint32_t and size_t, includes the SPIR-V array
    // (strip away directory names), and adds the sizeof.
    std::string c_file_name = args.preprocessed_file_name + ".c";

    size_t split_idx = 0;
    for (size_t i = 0; i < inc_file_name.size(); ++i) {
        if (inc_file_name[i] == '/') split_idx = i+1;
        #ifdef _WIN32
            if (inc_file_name[i] == '\\') split_idx = i+1;
        #endif
    }

    std::ofstream stream(c_file_name);
    stream << "#include <stddef.h>\n"
           << "#include <stdint.h>\n"
           << "#include \"" << &inc_file_name[split_idx] << "\"\n"
           << "size_t sizeof_" << variable_name
           << " = sizeof("<< variable_name << ");\n";
    if (!stream.is_open() || stream.bad() || stream.fail()) {
        const char* reason = strerror(errno);
        printf(
            "Could not write \"%s\": %s.\n",
            c_file_name.c_str(), reason);
        fprintf(
            stderr,
            "Could not write \"%s\": %s.\n",
            c_file_name.c_str(), reason);
        exit(1);
    }

    // Parent: wait for glsl compiler.
    int wstatus;
    auto waitpid_code = waitpid(pid, &wstatus, 0);
    if (waitpid_code < 0) {
        const char* msg = strerror(errno);
        printf("waitpid failed: %s.\n", msg);
        fprintf(stderr, "waitpid failed: %s.\n", msg);
        exit(1);
    }
    if (wstatus != 0) {
        printf("GLSL compiler failed: status %i.\n", wstatus);
        fprintf(stderr, "GLSL compiler failed: status %i.\n", wstatus);
        exit(WEXITSTATUS(wstatus) || 1);
    }
}

// Use execvp (replaces current process) to compile the preprocessed
// source file to the target file. Detect assembly or object file
// based on extension, and add -S or -c argument as appropriate.
//
// This is just one step for C/C++ files; for glsl, we have to use the
// glsl compiler to convert from glsl to C, then invoke the C
// compiler with the intermediate C file.
void exec_compile_to_target(Args& args)
{
    std::string c_input_filename = args.preprocessed_file_name;
    if (args.source_file_type == FileType::glsl) {
        c_input_filename += ".c";
        compile_glsl_to_c(args);
    }

    const char* cxx_arg = args.cxx.c_str();
    std::vector<const char*> argv;
    argv.push_back(cxx_arg);
    for (const auto& arg : args.cxxflags) {
        argv.push_back(arg.c_str());
    }
    argv.push_back(c_input_filename.c_str());

    // Just check the last character -- more thorough checking was
    // done during arg parsing.
    bool is_asm = (args.target_file_name.back() == 's');

    const char _s[] = "-S";
    const char _c[] = "-c";
    argv.push_back(is_asm ? _s : _c);

    const char _o[] = "-o";
    argv.push_back(_o);
    argv.push_back(args.target_file_name.c_str());

    const char* format = "%s";
    for (const char* arg : argv) {
        printf(format, arg);
        format = " %s";
    }
    printf("\n");

    argv.push_back(nullptr);
    fflush(stdout);
    fflush(stderr);
    execvp(cxx_arg, const_cast<char**>(argv.data()));

    const char* msg = strerror(errno);
    printf("Failed to execute compiler: %s\n", msg);
    fprintf(stderr, "Failed to execute compiler: %s\n", msg);
    exit(1);
}

int cckiss_main(Args& args)
{
    setbuf(stdout, nullptr);
    setbuf(stderr, nullptr);
    bool recompile = args.always_make || should_recompile_target_file(args);

    if (recompile) {
        preprocess_and_make_deps_file(args);

        std::string reason =
            args.always_make ? "-B passed to make" :
            args.changed_dependency_name.size() == 0 ? "not yet built" :
            ("change in " BLUE ) + args.changed_dependency_name + END;
        fprintf(stderr,
            "Compiling " BLUE "%s" END "; %s\n",
            args.target_file_name.c_str(),
            reason.c_str());

        exec_compile_to_target(args);
    }
    else {
        printf("#\"%s\": no dependency changes detected.\n",
            args.target_file_name.c_str());
        fprintf(stderr,
            "Up to date: " CYAN "%s\n" END,
            args.target_file_name.c_str());
    }

    return 0;
}

// Fake main used for rtags (rc -c -). Basically, print out a compile
// command that would have the same effect if cckiss were run for
// real, minus the intermediate preprocess stage (which confuses
// rtags).
int cckiss_rtags_main(Args& args)
{
    printf("%s", args.cxx.c_str());
    for (std::string& arg : args.cppflags) {
        printf(" %s", arg.c_str());
    }
    for (std::string& arg : args.cxxflags) {
        printf(" %s", arg.c_str());
    }
    bool is_asm = (args.target_file_name.back() == 's');
    printf(is_asm ? " -S" : " -c");

    printf(" %s -o %s\n",
        args.source_file_name.c_str(), args.target_file_name.c_str());
    return 0;
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    std::string delims =
        cxx_delim + ' ' + cxxflags_delim + ' ' + cppflags_delim
      + ' ' + glslc_delim + ' ' + glslargs_delim;

    Args args;

    const char* verbose_str = getenv("VERBOSE");
    args.verbose = verbose_str != nullptr && atoll(verbose_str) != 0;

    const char* makeargs_str = getenv("MAKEFLAGS");
    while (makeargs_str && *makeargs_str != '\0') {
        switch (*makeargs_str++) {
            case 'd': args.verbose = true; break;
            case 'B': args.always_make = true; break;
            case ' ': goto end_makeargs;
        }
    }
  end_makeargs:

    constexpr int cxxflags_mode = 0, cppflags_mode = 1, cxx_mode = 2,
                  glslc_mode = 3, glslargs_mode = 4;
    int mode = cxx_mode;

    if (argc >= 2) {
        args.target_file_name = argv[1];
    } else {
        printf("%s: First arg must be target file name.\n", argv[0]);
        fprintf(stderr, "%s: First arg must be target file name.\n", argv[0]);
        exit(1);
    }

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == cxxflags_delim) {
            mode = cxxflags_mode;
            continue;
        }
        else if (arg == cppflags_delim) {
            mode = cppflags_mode;
            continue;
        }
        else if (arg == cxx_delim) {
            mode = cxx_mode;
            continue;
        }
        else if (arg == glslc_delim) {
            mode = glslc_mode;
            continue;
        }
        else if (arg == glslargs_delim) {
            mode = glslargs_mode;
            continue;
        }
        else if (skip_prefix(arg.c_str(), ".cckiss.") != nullptr) {
            printf(
                "%s: Arg \"%s\" shall not start with '.cckiss.', "
                "unless it is a delimeter (%s).\n",
                argv[0], argv[i], delims.c_str());
            fprintf(
                stderr,
                "%s: Arg \"%s\" shall not start with '.cckiss.', "
                "unless it is a delimeter (%s).\n",
                argv[0], argv[i], delims.c_str());
            exit(1);
        }

        switch (mode) {
          default:
            assert(!"Unhandled mode in cckiss arg parsing");
          break; case cxxflags_mode:
            args.cxxflags.push_back(std::move(arg));
          break; case cppflags_mode:
            args.cppflags.push_back(std::move(arg));
          break; case cxx_mode:
            args.cxx = std::move(arg);
          break; case glslc_mode:
            args.glslc = std::move(arg);
          break; case glslargs_mode:
            args.glslargs.push_back(std::move(arg));
        }
    }

    // Initialize args struct from command line arguments.
    args.source_file_name = source_file_name_from_target(args.target_file_name);

    // This should be the one place where we determine the source file type
    // from file extension.
    const std::string& src_name = args.source_file_name;

    const auto sz = src_name.size();
    if (sz >= 5 && strcmp(&src_name[sz - 5], ".glsl") == 0) {
        args.source_file_type = FileType::glsl;
    }
    else if (sz >= 2 && strcmp(&src_name[sz - 2], ".c") == 0) {
        args.source_file_type = FileType::c;
    }
    else {
        args.source_file_type = FileType::cxx;
    }

    args.deps_file_name = deps_file_name_for_source(args.source_file_name);
    args.preprocessed_file_name = preprocessed_file_name_for_source(
        args.source_file_name,
        args.source_file_type);

    // Check for missing args.
    if (args.cxx == "") {
        printf(
            "%s: Blank CC or CXX arg (delim with %s).\n",
            argv[0], cxx_delim.c_str());
        fprintf(
            stderr,
            "%s: Blank CC or CXX arg (delim with %s).\n",
            argv[0], cxx_delim.c_str());
        exit(1);
    }
    if (args.glslc == "" && args.source_file_type == FileType::glsl) {
        printf(
            "%s: Blank GLSLC arg (delim with %s).\n",
            argv[0], glslc_delim.c_str());
        fprintf(
            stderr,
            "%s: Blank GLSLC arg (delim with %s).\n",
            argv[0], glslc_delim.c_str());
        exit(1);
    }

    if (getenv("CCKISS_RTAGS_PRINT")) {
        if (args.verbose) {
            fprintf(stderr, "Only printing compile command due "
                            "to CCKISS_RTAGS_PRINT env variable.\n");
        }
        return cckiss_rtags_main(args);
    }
    else {
        return cckiss_main(args);
    }
}
