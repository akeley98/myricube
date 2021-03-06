// Non-threadsafe shader loading & compilation helper. This is mostly
// what you would expect, except that there's this additional business
// of prepending some #defines (based on constexpr values in the C++
// side of things) so that shaders are automatically updated if
// parameters are tweaked.

#ifndef MYRICUBE_SHADERS_HH_
#define MYRICUBE_SHADERS_HH_

#include "myricube.hh"

#include <assert.h>
#include <errno.h>
#include <initializer_list>
#include <stdio.h>
#include <string>
#include <string.h>
#include <vector>

#include "glad/glad.h"
#include "MeshVoxelVertex.hh"
#include "PackedAABB.hh"

namespace myricube {

// Vertex shader input indices.
constexpr int packed_color_idx = 0;
constexpr int packed_vertex_idx = 1;
constexpr int packed_aabb_low_idx = 0;
constexpr int packed_aabb_high_idx = 1;

// Book-keeping for my plan to deliver arrays of voxel data to the GPU
// using shader storage buffer objects.

// glShaderStorageBlock argument 2 (storageBlockIndex).
constexpr int chunk_group_voxels_program_index = 0;

// glShaderStorageBlock argument 3 (storageBlockBinding).
constexpr int chunk_group_voxels_binding_index = 0;

// Data layout of said SSBO.
struct ChunkGroupVoxels
{
    uint32_t voxel_colors
        // chunk z/y/x in chunk group.
        [edge_chunks][edge_chunks][edge_chunks]
        // voxel z/y/x in chunk.
        [chunk_size][chunk_size][chunk_size];
};

// Image unit used by the staging compute shader, for the texture
// image being written to.
constexpr int staging_image_unit = 0;

// Binding point in the program for said image.
constexpr int staging_image_program_index = 0;

// Return a vector of #define lines and stuff (not newline terminated).
inline std::vector<std::string> get_preamble(std::string filename)
{
    return {
        "#version 450",
        "// FILE: " + filename,
        "#define CHUNK_SIZE " + std::to_string(chunk_size),
        "#define GROUP_SIZE " + std::to_string(group_size),
        // "#define BORDER_WIDTH " + std::to_string(border_width),
        "#define PACKED_VERTEX_IDX " + std::to_string(packed_vertex_idx),
        "#define PACKED_COLOR_IDX " + std::to_string(packed_color_idx),
        "#define PACKED_AABB_LOW_IDX " + std::to_string(packed_aabb_low_idx),
        "#define PACKED_AABB_HIGH_IDX " + std::to_string(packed_aabb_high_idx),
        "#define FOG_SCALAR 1.125",
        "#define RED_SHIFT " + std::to_string(red_shift),
        "#define GREEN_SHIFT " + std::to_string(green_shift),
        "#define BLUE_SHIFT " + std::to_string(blue_shift),
        uint32_t(visible_bit) == uint32_t(1) << 31 ?
            "#define VISIBLE_BIT int(1u << 31)" :
            "#define VISIBLE_BIT " + std::to_string(visible_bit),
        "#define X_SHIFT " + std::to_string(x_shift),
        "#define Y_SHIFT " + std::to_string(y_shift),
        "#define Z_SHIFT " + std::to_string(z_shift),
        "#define POS_X_FACE_BIT " + std::to_string(pos_x_face_bit),
        "#define NEG_X_FACE_BIT " + std::to_string(neg_x_face_bit),
        "#define POS_Y_FACE_BIT " + std::to_string(pos_y_face_bit),
        "#define NEG_Y_FACE_BIT " + std::to_string(neg_y_face_bit),
        "#define POS_Z_FACE_BIT " + std::to_string(pos_z_face_bit),
        "#define NEG_Z_FACE_BIT " + std::to_string(neg_z_face_bit),
        "#define CHUNK_GROUP_VOXELS_PROGRAM_INDEX " +
            std::to_string(chunk_group_voxels_program_index),
        "#define STAGING_IMAGE_PROGRAM_INDEX " +
            std::to_string(staging_image_program_index),
    };
}

// Compile and link a shader using the shader files in the list of
// `filename_count` strings. (They are found in the
// data_directory). Vertex shader files end with .vert", geometry
// shaders with ".geom", and fragment shaders with ".frag".
inline GLuint make_program(const char* const* filenames, size_t filename_count);

inline GLuint make_program(std::initializer_list<const char*> filenames)
{
    return make_program(&*filenames.begin(), filenames.size());
}

// Return the full source code (as a string) of the shader in the
// named file, prepending the stuff from get_preamble and
// a comment indicating its filename.
inline std::string read_shader_source(const char* raw_filename) noexcept
{
#ifdef MYRICUBE_WINDOWS
    filename_string actual_filename = expand_filename(raw_filename);
    FILE* file = _wfopen(actual_filename.c_str(), L"r");
    std::string filename;
    for (filename_char c : actual_filename) {
        filename.push_back(c <= 127 ? c : '?'); // IDGAF Windows is annoying
    }
#else
    std::string filename = expand_filename(raw_filename);
    FILE* file = fopen(filename.c_str(), "r");
#endif

    if (file == nullptr) {
        panic("Could not open " + filename + " (" + std::to_string(errno)
            + ") " + strerror(errno));
    }
    errno = 0;

    // I know there are better ways than just +='ing the whole file,
    // but it doesn't really matter at the moment. Shaders are not all
    // that big.
    std::string result;
    auto defines = get_preamble(filename);

    // Replace blank lines in the start of the source code with the
    // #defines. This prevents the line numbers from being off.
    bool warn = true;
    for (size_t i = 0; i < defines.size(); ++i) {
        int c = getc(file);
        if (c != '\n') {
            if (warn) {
                fprintf(stderr, "Line numbers in %s might be off "
                    "due to prepended #defines.\n", filename.c_str());
                warn = false;
            }
            ungetc(c, file);
        }
        result += defines[i] + "\n";
    }

    int c;
    while ((c = getc(file)) != EOF) {
        result += char(c);
    }

    if (errno != 0) {
        std::string message = "Warning: Unexpected EOF: " + filename + " (" +
            std::to_string(errno) + ") " + strerror(errno);
        fprintf(stderr, "%s\n", message.c_str());
    }
    return result;
}

// Print to stderr the source code with line numbers.
inline void print_source(const std::string& source)
{
    const char format[] = "%4i ";

    int line_number = 1;
    fprintf(stderr, format, 1);
    for (char c : source)
    {
        putc(c, stderr);
        if (c == '\n') {
            fprintf(stderr, format, ++line_number);
        }
    }
    putc('\n', stderr);
}

inline bool is_vertex_shader_filename(const char* filename) {
    auto len = strlen(filename);
    return len >= 5 and strcmp(".vert", &filename[len-5]) == 0;
}

inline bool is_geometry_shader_filename(const char* filename) {
    auto len = strlen(filename);
    return len >= 5 and strcmp(".geom", &filename[len-5]) == 0;
}

inline bool is_fragment_shader_filename(const char* filename) {
    auto len = strlen(filename);
    return len >= 5 and strcmp(".frag", &filename[len-5]) == 0;
}

inline bool is_compute_shader_filename(const char* filename) {
    auto len = strlen(filename);
    return len >= 5 and strcmp(".comp", &filename[len-5]) == 0;
}

inline GLuint make_program(const char* const* filenames, size_t filename_count)
{
    // Outputs for compiler messages.
    std::vector<char> log(100000);
    GLint okay = 0;
    GLsizei length = 0;

    // Create vectors of each shader type.
    std::vector<GLuint> vertex_shaders, geometry_shaders, fragment_shaders,
                        compute_shaders;

    // Load each shader file, compile, and place into the correct vector.
    for (size_t i = 0; i < filename_count; ++i) {
        const char* f = filenames[i];
        fprintf(stderr, "Loading shader %s...\n", f);

        std::string source = read_shader_source(f);
        const char* source_c_str = source.c_str();
        GLuint id = 0;

        if (is_vertex_shader_filename(f)) {
            id = glCreateShader(GL_VERTEX_SHADER);
            vertex_shaders.push_back(id);
        }
        else if (is_geometry_shader_filename(f)) {
            id = glCreateShader(GL_GEOMETRY_SHADER);
            geometry_shaders.push_back(id);
            // GS support not tested after rewrite!!!
        }
        else if (is_fragment_shader_filename(f)) {
            id = glCreateShader(GL_FRAGMENT_SHADER);
            fragment_shaders.push_back(id);
        }
        else if (is_compute_shader_filename(f)) {
            id = glCreateShader(GL_COMPUTE_SHADER);
            compute_shaders.push_back(id);
        }
        else {
            panic(f + std::string(" should end in .frag or .geom or .vert or .comp"));
        }
        glShaderSource(id, 1, &source_c_str, nullptr);

        PANIC_IF_GL_ERROR;

        glCompileShader(id);
        glGetShaderiv(id, GL_COMPILE_STATUS, &okay);
        glGetShaderInfoLog(id, log.size()-1, &length, log.data());

        if (!okay or length > 0) {
            print_source(source);
            fputs(log.data(), stderr);
        }
        if (!okay) {
            panic("Shader compile error: " + std::string(f));
        }
    }

    // Vertex and fragment shaders are mandatory unless there is a
    // compute shader.
    if (compute_shaders.size() == 0) {
        if (vertex_shaders.size() == 0) {
            panic("No vertex shaders.");
        }
        if (fragment_shaders.size() == 0) {
            panic("No fragment shaders.");
        }
    }

    // Create a program object and put each shader in it.
    GLuint program_id = glCreateProgram();
    for (auto id : vertex_shaders) {
        glAttachShader(program_id, id);
    }
    for (auto id : geometry_shaders) {
        glAttachShader(program_id, id);
    }
    for (auto id : fragment_shaders) {
        glAttachShader(program_id, id);
    }
    for (auto id : compute_shaders) {
        glAttachShader(program_id, id);
    }
    PANIC_IF_GL_ERROR;

    // Link and return the program.
    glLinkProgram(program_id);
    glGetProgramiv(program_id, GL_LINK_STATUS, &okay);
    glGetProgramInfoLog(program_id, log.size()-1, &length, log.data());

    if (!okay or length > 0) {
        fprintf(stderr, "%s\n", log.data());
    }
    if (!okay) {
        panic("Shader link error");
    }
    PANIC_IF_GL_ERROR;
    fprintf(stderr, "...program id=%i\n", int(program_id));
    return program_id;
}

} // end namespace

#endif /* !MYRICUBE_SHADERS_HH_ */
