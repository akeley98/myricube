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

namespace myricube {

constexpr float border_width = 0.1f;

// Vertex shader input indices.
constexpr int packed_color_idx = 0;
constexpr int packed_vertex_idx = 1;

constexpr int packed_aabb_low_idx = 0;
constexpr int packed_aabb_high_idx = 1;
constexpr int unit_box_vertex_idx = 2;
constexpr int unit_box_normal_idx = 3;

// bit assignments for packed verts.
constexpr int x_shift = 0;
constexpr int y_shift = 8;
constexpr int z_shift = 16;

constexpr int pos_x_face_bit = (1 << 24);
constexpr int neg_x_face_bit = (1 << 25);
constexpr int pos_y_face_bit = (1 << 26);
constexpr int neg_y_face_bit = (1 << 27);
constexpr int pos_z_face_bit = (1 << 28);
constexpr int neg_z_face_bit = (1 << 29);

// Return a vector of #define lines and stuff (not newline terminated).
inline std::vector<std::string> get_preamble(std::string filename)
{
    return {
        "#version 450",
        "// FILE: " + filename,
        "#define CHUNK_SIZE " + std::to_string(chunk_size),
        "#define GROUP_SIZE " + std::to_string(group_size),
        "#define BORDER_WIDTH " + std::to_string(border_width),
        "#define PACKED_VERTEX_IDX " + std::to_string(packed_vertex_idx),
        "#define PACKED_COLOR_IDX " + std::to_string(packed_color_idx),
        "#define UNIT_BOX_VERTEX_IDX " + std::to_string(unit_box_vertex_idx),
        "#define UNIT_BOX_NORMAL_IDX " + std::to_string(unit_box_normal_idx),
        "#define PACKED_AABB_LOW_IDX " + std::to_string(packed_aabb_low_idx),
        "#define PACKED_AABB_HIGH_IDX " + std::to_string(packed_aabb_high_idx),
        "#define FOG_SCALAR 1.125",
        "#define X_SHIFT " + std::to_string(x_shift),
        "#define Y_SHIFT " + std::to_string(y_shift),
        "#define Z_SHIFT " + std::to_string(z_shift),
        "#define POS_X_FACE_BIT " + std::to_string(pos_x_face_bit),
        "#define NEG_X_FACE_BIT " + std::to_string(neg_x_face_bit),
        "#define POS_Y_FACE_BIT " + std::to_string(pos_y_face_bit),
        "#define NEG_Y_FACE_BIT " + std::to_string(neg_y_face_bit),
        "#define POS_Z_FACE_BIT " + std::to_string(pos_z_face_bit),
        "#define NEG_Z_FACE_BIT " + std::to_string(neg_z_face_bit),
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
    std::string filename = expand_filename(raw_filename);
    FILE* file = fopen(filename.c_str(), "r");
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

inline GLuint make_program(const char* const* filenames, size_t filename_count)
{
    // Output memory for compiler messages.
    std::vector<char> log(100000);

    // Store vertex, geometry, and fragment shader sources
    // separately. Use std::string to store shader source code (loaded
    // from files), but we also need a parallel const char* array for
    // the OpenGL C API.
    std::vector<std::string> vs_string_array;
    std::vector<std::string> gs_string_array;
    std::vector<std::string> fs_string_array;
    std::vector<const char*> vs_c_str_array;
    std::vector<const char*> gs_c_str_array;
    std::vector<const char*> fs_c_str_array;

    // Load each shader file and shunt the source code to the correct
    // array depending on its file extension.
    for (size_t i = 0; i < filename_count; ++i) {
        const char* f = filenames[i];

        if (is_vertex_shader_filename(f)) {
            vs_string_array.push_back(read_shader_source(f));
            vs_c_str_array.push_back(vs_string_array.back().data());
        }
        else if (is_geometry_shader_filename(f)) {
            gs_string_array.push_back(read_shader_source(f));
            gs_c_str_array.push_back(gs_string_array.back().data());
        }
        else if (is_fragment_shader_filename(f)) {
            fs_string_array.push_back(read_shader_source(f));
            fs_c_str_array.push_back(fs_string_array.back().data());
        }
        else {
            panic(f + std::string(" should end in .frag or .geom or .vert"));
        }
    }

    // Make OpenGL shader objects.
    GLuint program_id = glCreateProgram();
    GLuint vs_id = glCreateShader(GL_VERTEX_SHADER);
    GLuint fs_id = glCreateShader(GL_FRAGMENT_SHADER);
    GLuint gs_id = 0;

    // Mandatory vertex shader.
    if (vs_c_str_array.size() != 0) {
        glShaderSource(vs_id, vs_c_str_array.size(),
                       vs_c_str_array.data(), nullptr);
    }
    else {
        panic("No vertex shaders (.vert)");
    }

    // Mandatory fragment shader.
    if (fs_c_str_array.size() != 0) {
        glShaderSource(fs_id, fs_c_str_array.size(),
                       fs_c_str_array.data(), nullptr);
    }
    else {
        panic("No fragment shaders (.frag)");
    }

    // Optional geometry shader.
    if (gs_c_str_array.size() != 0) {
        gs_id = glCreateShader(GL_GEOMETRY_SHADER);
        glShaderSource(gs_id, gs_c_str_array.size(),
                       gs_c_str_array.data(), nullptr);
    }
    /// else {
    //     fprintf(stderr, "Note: no geometry shaders (.geom).\n");
    // }
    PANIC_IF_GL_ERROR;

    // Compile each shader.
    GLint okay = 0;
    GLsizei length = 0;
    const GLuint shader_id_array[3] = { vs_id, gs_id, fs_id };

    for (auto id : shader_id_array) {
        // Skip geometry shader if there is none.
        if (id == 0) continue;

        glCompileShader(id);
        glGetShaderiv(id, GL_COMPILE_STATUS, &okay);
        if (okay) {
            glAttachShader(program_id, id);
        }
        glGetShaderInfoLog(id, log.size()-1, &length, log.data());

        if (!okay or length > 0) {
            auto& sources = id == vs_id ? vs_string_array
                          : id == fs_id ? fs_string_array : gs_string_array;
            for (const std::string& source : sources) {
                print_source(source);
            }
            fputs(log.data(), stderr);
        }

        if (!okay) panic("Shader compile error.");
    }

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
    return program_id;
}

} // end namespace

#endif /* !MYRICUBE_SHADERS_HH_ */
