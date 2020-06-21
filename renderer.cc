// Implementation of the hybrid voxel renderer. I don't see how I can
// make much of this exception safe, so I wrap everything in a
// noexcept at the end. As usual with OpenGL, none of this is
// thread-safe either.

#include "myricube.hh"

#include <stdio.h>
#include <typeinfo>
#include <utility>
#include <vector>

#include "camera.hh"
#include "chunk.hh"
#include "glad/glad.h"
#include "renderer.hh"
#include "SDL2/SDL.h"

namespace myricube {

static constexpr int edge_chunks = group_size / chunk_size;
bool chunk_debug = false;

// TODO: Maybe get some better code for compiling shaders? Load it
// from a file?
static GLuint make_program(const char* vs_code, const char* fs_code)
{
    static GLchar log[1024];
    PANIC_IF_GL_ERROR;
    GLuint program_id = glCreateProgram();
    GLuint vs_id = glCreateShader(GL_VERTEX_SHADER);
    GLuint fs_id = glCreateShader(GL_FRAGMENT_SHADER);

    const GLchar* string_array[1];
    string_array[0] = (GLchar*)vs_code;
    glShaderSource(vs_id, 1, string_array, nullptr);
    string_array[0] = (GLchar*)fs_code;
    glShaderSource(fs_id, 1, string_array, nullptr);

    glCompileShader(vs_id);
    glCompileShader(fs_id);

    PANIC_IF_GL_ERROR;

    GLint okay = 0;
    GLsizei length = 0;
    const GLuint shader_id_array[2] = { vs_id, fs_id };
    for (auto id : shader_id_array) {
        glGetShaderiv(id, GL_COMPILE_STATUS, &okay);
        if (okay) {
            glAttachShader(program_id, id);
        } else {
            glGetShaderInfoLog(id, sizeof log, &length, log);
            fprintf(stderr, "%s\n", id == vs_id ? vs_code : fs_code);
            panic("Shader compilation error", log);
        }
    }

    glLinkProgram(program_id);
    glGetProgramiv(program_id, GL_LINK_STATUS, &okay);
    if (!okay) {
        glGetProgramInfoLog(program_id, sizeof log, &length, log);
        panic("Shader link error", log);
    }

    PANIC_IF_GL_ERROR;
    return program_id;
}

// Single vertex of the VBO for a chunk's mesh.
struct MeshVertex
{
    uint32_t packed_vertex;
    uint32_t packed_color;

    // Pack this vertex given the vertex's residue coordinates
    // (i.e. local coordinate in coordinate system of a chunk group)
    // and 24-bit color.
    //
    // This constructor defines the packed layout.
    MeshVertex(uint8_t x, uint8_t y, uint8_t z,
               uint8_t red, uint8_t green, uint8_t blue)
    {
        static_assert(group_size <= 255,
                      "group too big for 8-bit unsigned coordinates.");
        packed_vertex = x | uint32_t(y) << 8 | uint32_t(z) << 16;
        packed_color = blue | uint32_t(green) << 8 | uint32_t(red) << 16;
    }
};

// CPU-side copy of the mesh for a single chunk.
struct ChunkMesh
{
    // CPU vertex list.
    std::vector<MeshVertex> verts;

    // Offset of first byte in the VBO where the GPU copy of the
    // vertex list is stored.
    GLsizeiptr vbo_byte_offset;

    // Number of bytes reserved in the VBO for this chunk.
    GLsizeiptr vbo_bytes_reserved;
};

// Handle for GPU resources for the mesh of one chunk group.  The
// CPU-side state of this MeshEntry should always accurately describe
// the state of the copy of this data on the GPU and GL. (i.e. if you
// modify the data within, it is your responsibility to update the GPU
// state at the same time).
struct MeshEntry
{
    // World ID and group coordinate of the chunk group this mesh is
    // for. This can be used to tell if this entry is correct or needs
    // to be replaced in MeshStore.
    uint64_t world_id = uint64_t(-1);
    glm::ivec3 group_coord;

    // mesh_array[z][y][x] is the mesh for ChunkGroup::chunk_array[z][y][x].
    ChunkMesh mesh_array[edge_chunks][edge_chunks][edge_chunks];

    // OpenGL name for the VBO used to store the meshes of this chunk group.
    // 0 when not yet allocated.
    GLuint vbo_name = 0;

    // Size in bytes of the vbo's data store on the GPU.
    GLsizeiptr vbo_bytes = 0;

    // Bytes used for the VBO (i.e. the byte offset where new data can be
    // written to the vbo).
    GLsizeiptr bytes_used = 0;

    MeshEntry() = default;

    // Get some RAII going.
    ~MeshEntry()
    {
        glDeleteBuffers(1, &vbo_name);
        vbo_name = 0;
        vbo_bytes = 0;
        bytes_used = 0;
        world_id = uint64_t(-1);
    }

    // Disable moves, but, see swap below. (Harkens back to C++98 hacks :P )
    MeshEntry(MeshEntry&&) = delete;
};

void swap(MeshEntry& left, MeshEntry& right)
{
    using std::swap;
    swap(left.world_id, right.world_id);
    swap(left.group_coord, right.group_coord);
    swap(left.mesh_array, right.mesh_array);
    swap(left.vbo_name, right.vbo_name);
    swap(left.vbo_bytes, right.vbo_bytes);
    swap(left.bytes_used, right.bytes_used);
}

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
        packed_low = uint32_t(x) | uint32_t(y) << 8 | uint32_t(z) << 16;
        x = aabb_high.x;
        y = aabb_high.y;
        z = aabb_high.z;
        packed_high = uint32_t(x) | uint32_t(y) << 8 | uint32_t(z) << 16;
    }
};

// All voxels within a chunk group share a single 3D texture on the
// GPU, as well as one VBO used to store the array of AABB for the
// chunks within the group. This is the handle for the GPU data needed
// to raycast one chunk group.
struct RaycastEntry
{
    // World ID and group coordinate of the chunk group this texture
    // and AABB array is for. This can be used to tell if this entry
    // is correct or needs to be replaced in RaycastStore.
    uint64_t world_id = uint64_t(-1);
    glm::ivec3 group_coord;

    // aabb_array[z][y][x] is the minimal AABB containing the visible
    // voxels of ChunkGroup::chunk_array[z][y][x].
    //
    // If you change this, be aware that sizeof is used on this to
    // know the amount of bytes to allocate for the GPU's VBO.
    PackedAABB aabb_array[edge_chunks][edge_chunks][edge_chunks];

    // OpenGL name of the VBO storing the aabb_array. 0 when not yet
    // allocated.
    GLuint vbo_name = 0;

    // OpenGL name for the 3D texture representing this group of
    // voxels on the GPU. 0 when not yet allocated.
    GLuint texture_name = 0;

    RaycastEntry() = default;

    ~RaycastEntry()
    {
        glDeleteBuffers(1, &vbo_name);
        glDeleteTextures(1, &texture_name);
        vbo_name = 0;
        texture_name = 0;
        world_id = uint64_t(-1);
    }

    RaycastEntry(RaycastEntry&&) = delete;
};

void swap(RaycastEntry& left, RaycastEntry& right)
{
    using std::swap;
    swap(left.world_id, right.world_id);
    swap(left.group_coord, right.group_coord);
    swap(left.aabb_array, right.aabb_array);
    swap(left.vbo_name, right.vbo_name);
    swap(left.texture_name, right.texture_name);
}

// Cache of MeshEntry/RaycastEntry objects. The idea is to store entry
// structs for chunk groups near the camera. They are stored in a
// wraparound fashion (i.e.  using modular arithmetic on group
// coordinates) so that as the camera moves, the new groups being
// rendered overwrites the old groups no longer being rendered.
template <class Entry, int N>
class BaseStore
{
    // N should be made tunable later depending on
    // Camera::raycast_threshold or Camera::far_plane.

    Entry entry_array[N][N][N];

  public:
    // Return a pointer to the location that the Entry for the chunk
    // group with the given group coordinate would be, IF it exists.
    // (i.e. that location may be blank or occupied by an entry for
    // another group).
    Entry* cached_location(glm::ivec3 group_coord)
    {
        uint32_t x = (uint32_t(group_coord.x) ^ 0x1000'0000) % N;
        uint32_t y = (uint32_t(group_coord.y) ^ 0x1000'0000) % N;
        uint32_t z = (uint32_t(group_coord.z) ^ 0x1000'0000) % N;
        return &entry_array[z][y][x];
    }

    // Return a valid, updated Entry for the given chunk group
    // extracted from the given world. TODO document templates needed
    // for this to work.
    template <typename EntryFiller>
    Entry* update(PositionedChunkGroup& pcg, VoxelWorld& world)
    {
        glm::ivec3 gc = group_coord(pcg);
        Entry* entry = cached_location(gc);

        // If the Entry in the array already corresponds to the given
        // chunk group; just update it (deal with dirty chunks) and
        // return.
        if (entry->world_id == world.id()
        and entry->group_coord == group_coord(pcg)) {
            EntryFiller::update(pcg, world, entry);
            return entry;
        }

        // Otherwise, we need to evict the current entry and repopulate it.
        if (long(entry->world_id) >= 0) {
            fprintf(stderr, "Evicting %s entry for chunk group "
                            "(%i %i %i) of world %li\n",
                            typeid(Entry).name(),
                            int(entry->group_coord.x),
                            int(entry->group_coord.y),
                            int(entry->group_coord.z),
                            long(entry->world_id));

        }
        EntryFiller::replace(pcg, world, entry);
        return entry;
    }
};

class MeshStore : public BaseStore<MeshEntry, 5> { };
class RaycastStore : public BaseStore<RaycastEntry, 7> { };

#define BORDER_WIDTH_STR "0.1"

const char mesh_vs_source[] =
"#version 330\n"
"layout(location=0) in int packed_vertex;\n"
"layout(location=1) in int packed_color;\n"
"out vec3 color;\n"
"out vec3 model_space_position_;\n"
"uniform mat4 mvp_matrix;\n"
"void main() {\n"
    "float x = float(packed_vertex & 255);\n"
    "float y = float((packed_vertex >> 8) & 255);\n"
    "float z = float((packed_vertex >> 16) & 255);\n"
    "vec4 model_space_position = vec4(x, y, z, 1);\n"
    "gl_Position = mvp_matrix * model_space_position;\n"
    "float red   = ((packed_color >> 16) & 255) * (1./255.);\n"
    "float green = ((packed_color >> 8) & 255) * (1./255.);\n"
    "float blue  = (packed_color & 255) * (1./255.);\n"
    "color = vec3(red, green, blue);\n"
    "model_space_position_ = model_space_position.xyz;\n"
"}\n";

constexpr int packed_vertex_idx = 0;
constexpr int packed_color_idx = 1;

const char mesh_fs_source[] =
"#version 330\n"
"in vec3 color;\n"
"in vec3 model_space_position_;\n"
"out vec4 out_color;\n"
"void main() {\n"
    "const float d = " BORDER_WIDTH_STR ";\n"
    "float x = model_space_position_.x;\n"
    "int x_border = (x - floor(x + d) < d) ? 1 : 0;\n"
    "float y = model_space_position_.y;\n"
    "int y_border = (y - floor(y + d) < d) ? 1 : 0;\n"
    "float z = model_space_position_.z;\n"
    "int z_border = (z - floor(z + d) < d) ? 1 : 0;\n"
    "float scale = (x_border + y_border + z_border >= 2) ? 0.5 : 1.0;\n"
    "out_color = vec4(color * scale, 1);\n"
"}\n";

// Really need to improve this...
#define GROUP_SIZE_STR "64"
constexpr int unit_box_vertex_idx = 0;
constexpr int packed_aabb_low_idx = 1;
constexpr int packed_aabb_high_idx = 2;

static const char raycast_vs_source[] =
"#version 330\n"
"layout(location=0) in vec4 unit_box_vertex;\n"
"layout(location=1) in int packed_aabb_low;\n"
"layout(location=2) in int packed_aabb_high;\n"
"out vec3 residue_coord;\n"
"out float border_fade;\n"
"flat out ivec3 aabb_low;\n"
"flat out ivec3 aabb_high;\n"
"uniform mat4 mvp_matrix;\n"
"uniform vec3 eye_relative_group_origin;\n"
"void main() {\n"
    "int low_x = packed_aabb_low & 255;\n"
    "int low_y = (packed_aabb_low >> 8) & 255;\n"
    "int low_z = (packed_aabb_low >> 16) & 255;\n"
    "int high_x = packed_aabb_high & 255;\n"
    "int high_y = (packed_aabb_high >> 8) & 255;\n"
    "int high_z = (packed_aabb_high >> 16) & 255;\n"
    "vec3 f_aabb_low = vec3(low_x, low_y, low_z);\n"
    "vec3 sz = vec3(high_x, high_y, high_z) - f_aabb_low;\n"
    "vec4 model_space_pos = vec4(unit_box_vertex.xyz * sz + f_aabb_low, 1);\n"
    "gl_Position = mvp_matrix * model_space_pos;\n"
    "vec3 disp = model_space_pos.xyz - eye_relative_group_origin;\n"
    "float distance = sqrt(dot(disp, disp));\n"
    "border_fade = clamp(sqrt(distance) * 0.042, 0.5, 1.0);\n"
    "aabb_low = ivec3(low_x, low_y, low_z);\n"
    "aabb_high = ivec3(high_x, high_y, high_z);\n"
    "residue_coord = model_space_pos.xyz;\n"
"}\n";

static const char raycast_fs_source[] =
"#version 330\n"
"in vec3 residue_coord;\n"
"in float border_fade;\n"
"flat in ivec3 aabb_low;\n"
"flat in ivec3 aabb_high;\n"
// TODO: Add bias (based on normal vector) to deal with floor/ceil
// rounding errors. Also hide this mess in a glsl file somewhere.
"uniform vec3 eye_relative_group_origin;\n"
"uniform sampler3D chunk_blocks;\n"
"uniform bool chunk_debug;\n"
"out vec4 color;\n"
"void main() {\n"
"if (!chunk_debug) {\n"
    "const float d = " BORDER_WIDTH_STR ";\n"
    "float x0 = eye_relative_group_origin.x;\n"
    "float y0 = eye_relative_group_origin.y;\n"
    "float z0 = eye_relative_group_origin.z;\n"
    "vec3 slope = vec3(residue_coord) - eye_relative_group_origin;\n"
    "float xm = slope.x;\n"
    "float ym = slope.y;\n"
    "float zm = slope.z;\n"
    "float rcp = 1.0/" GROUP_SIZE_STR ".0;\n"
    "float best_t = 1.0 / 0.0;\n"
    "vec4 best_color = vec4(0,0,0,0);\n"
    "vec3 best_coord = vec3(0,0,0);\n"
    "int iter = 0;\n"
    "int x_init = int(xm > 0 ? ceil(residue_coord.x) \n"
                            ": floor(residue_coord.x));\n"
    "int x_end = xm > 0 ? aabb_high.x : aabb_low.x;\n"
    // "int x_end = xm > 0 ? 64 : 0;\n"
    "int x_step = xm > 0 ? 1 : -1;\n"
    "float x_fudge = xm > 0 ? .25 : -.25;\n"
    "for (int x = x_init; x != x_end; x += x_step) {\n"
        "if (iter++ >= 255) { color = vec4(1,0,1,1); return; }\n"
        "float t = (x - x0) / xm;\n"
        "float y = y0 + ym * t;\n"
        "float z = z0 + zm * t;\n"
        "if (y < aabb_low.y || y > aabb_high.y) break;\n"
        "if (z < aabb_low.z || z > aabb_high.z) break;\n"
        "vec3 texcoord = vec3(x + x_fudge, y, z) * rcp;\n"
        "vec4 lookup_color = texture(chunk_blocks, texcoord);\n"
        "if (lookup_color.a > 0 && t > 0) {\n"
            "if (best_t > t) {\n"
                "best_t = t;\n"
                "best_color = lookup_color;\n"
                "best_coord = vec3(x,y,z);\n"
                "if (y - floor(y + d) < d || z - floor(z + d) < d) {\n"
                    "best_color.rgb *= border_fade;\n"
                "}\n"
            "}\n"
            "break;\n"
        "}\n"
    "}\n"
    "int y_init = int(ym > 0 ? ceil(residue_coord.y) \n"
                            ": floor(residue_coord.y));\n"
    "int y_end = ym > 0 ? aabb_high.y : aabb_low.y;\n"
    // "int y_end = ym > 0 ? 64 : 0;\n"
    "int y_step = ym > 0 ? 1 : -1;\n"
    "float y_fudge = ym > 0 ? .25 : -.25;\n"
    "for (int y = y_init; y != y_end; y += y_step) {\n"
        "if (iter++ >= 255) { color = vec4(1,0,1,1); return; }\n"
        "float t = (y - y0) / ym;\n"
        "float x = x0 + xm * t;\n"
        "float z = z0 + zm * t;\n"
        "if (x < aabb_low.x || x > aabb_high.x) break;\n"
        "if (z < aabb_low.z || z > aabb_high.z) break;\n"
        "vec3 texcoord = vec3(x, y + y_fudge, z) * rcp;\n"
        "vec4 lookup_color = texture(chunk_blocks, texcoord);\n"
        "if (lookup_color.a > 0 && t > 0) {\n"
            "if (best_t > t) {\n"
                "best_t = t;\n"
                "best_color = lookup_color;\n"
                "best_coord = vec3(x,y,z);\n"
                "if (x - floor(x + d) < d || z - floor(z + d) < d) {\n"
                    "best_color.rgb *= border_fade;\n"
                "}\n"
            "}\n"
            "break;\n"
        "}\n"
    "}\n"
    "int z_init = int(zm > 0 ? ceil(residue_coord.z) \n"
                            ": floor(residue_coord.z));\n"
    "int z_end = zm > 0 ? aabb_high.z : aabb_low.z;\n"
    // "int z_end = zm > 0 ? 64 : 0;\n"
    "int z_step = zm > 0 ? 1 : -1;\n"
    "float z_fudge = zm > 0 ? .25 : -.25;\n"
    "for (int z = z_init; z != z_end; z += z_step) {\n"
        "if (iter++ >= 255) { color = vec4(1,0,1,1); return; }\n"
        "float t = (z - z0) / zm;\n"
        "float x = x0 + xm * t;\n"
        "float y = y0 + ym * t;\n"
        "if (x < aabb_low.x || x > aabb_high.x) break;\n"
        "if (y < aabb_low.y || y > aabb_high.y) break;\n"
        "vec3 texcoord = vec3(x, y, z + z_fudge) * rcp;\n"
        "vec4 lookup_color = texture(chunk_blocks, texcoord);\n"
        "if (lookup_color.a > 0 && t > 0) {\n"
            "if (best_t > t) {\n"
                "best_t = t;\n"
                "best_color = lookup_color;\n"
                "best_coord = vec3(x,y,z);\n"
                "if (x - floor(x + d) < d || y - floor(y + d) < d) {\n"
                    "best_color.rgb *= border_fade;\n"
                "}\n"
            "}\n"
            "break;\n"
        "}\n"
    "}\n"
    "if (best_color.a == 0) discard;\n"
    "color = best_color;\n"
    //"vec4 v = vp_matrix * vec4(best_coord + chunk_offset, 1);\n"
    //"gl_FragDepth = .3;\n"
"} else {\n"
    "int x_floor = int(floor(residue_coord.x));\n"
    "int y_floor = int(floor(residue_coord.y));\n"
    "int z_floor = int(floor(residue_coord.z));\n"
    "color = vec4(x_floor & 1, y_floor & 1, z_floor & 1, 1);\n"
"}\n"
"}\n";

static const float unit_box_vertices[32] =
{
    0, 1, 1, 1,
    0, 0, 1, 1,
    1, 0, 1, 1,
    1, 1, 1, 1,
    0, 1, 0, 1,
    0, 0, 0, 1,
    1, 0, 0, 1,
    1, 1, 0, 1,
};

static const GLushort unit_box_elements[36] = {
    6, 7, 2, 7, 3, 2,
    4, 5, 0, 5, 1, 0,
    0, 3, 4, 3, 7, 4,
    6, 2, 5, 2, 1, 5,
    2, 3, 1, 3, 0, 1,
    6, 5, 7, 5, 4, 7,
};

// Renderer class. Instantiate it with the camera and world to render
// and use it once.
class Renderer
{
    Camera& camera;
    VoxelWorld& world;
  public:
    Renderer(VoxelWorld& world_, Camera& camera_) :
        camera(camera_), world(world_) { }

    // Given a chunk, make the mesh needed to render it. I also need
    // the residue coordinate of the lower-left of the chunk (i.e.
    // the position of the chunk relative to the chunk group it is
    // in) because the mesh is defined to be in residue coordinates.
    //
    // This is just ugly.
    //
    // TODO: consider properly handling edge cases with other chunks?
    // It's less easy than it seems, as the mesh can change not only
    // due to changes in this chunk, but also changes in neighbor
    // chunks that reveal previously hidden voxel faces. This may
    // not be worth it, so for now, I don't include a VoxelWorld
    // argument needed to make this happen.
    static std::vector<MeshVertex> make_mesh_verts(const Chunk& chunk,
                                                   glm::ivec3 chunk_residue)
    {
        // Look up whether the voxel at the given coordinate
        // (relative to the lower-left of this chunk) is visible.
        // Act as if voxels outside the chunk are always invisible.
        auto visible_block = [&] (glm::ivec3 coord) -> bool
        {
            if (coord.x < 0 or coord.x >= chunk_size
             or coord.y < 0 or coord.y >= chunk_size
             or coord.z < 0 or coord.z >= chunk_size) return false;
            // Note: the masking in Chunk::operator () won't mess up
            // this coord. (I'm kind of violating my own comment in
            // Chunk::operator() because coord won't actually be in
            // the chunk unless that chunk is at (0,0,0)).
            return chunk(coord).visible;
        };

        MeshVertex reference(
                chunk_residue.x,
                chunk_residue.y,
                chunk_residue.z,
                0, 0, 0);
        std::vector<MeshVertex> verts;

        auto verts_add_block = [&] (glm::ivec3 coord)
        {
            Voxel v = chunk(coord);
            if (!v.visible) return;

            // printf("%i %i %i\n", int(coord.x), int(coord.y), int(coord.z));

            auto r = v.red;
            auto g = v.green;
            auto b = v.blue;
            uint8_t x = uint8_t(coord.x) + chunk_residue.x;
            uint8_t y = uint8_t(coord.y) + chunk_residue.y;
            uint8_t z = uint8_t(coord.z) + chunk_residue.z;
            uint8_t x1 = x+1;
            uint8_t y1 = y+1;
            uint8_t z1 = z+1;

            if (!visible_block(coord + glm::ivec3(-1, 0, 0))) {
                verts.emplace_back(x,  y1, z1, r, g, b);
                verts.emplace_back(x,  y1, z,  r, g, b);
                verts.emplace_back(x,  y,  z1, r, g, b);
                verts.emplace_back(x,  y1, z,  r, g, b);
                verts.emplace_back(x,  y,  z,  r, g, b);
                verts.emplace_back(x,  y,  z1, r, g, b);
            }
            if (!visible_block(coord + glm::ivec3(1, 0, 0))) {
                verts.emplace_back(x1, y,  z,  r, g, b);
                verts.emplace_back(x1, y1, z,  r, g, b);
                verts.emplace_back(x1, y,  z1, r, g, b);
                verts.emplace_back(x1, y1, z,  r, g, b);
                verts.emplace_back(x1, y1, z1, r, g, b);
                verts.emplace_back(x1, y,  z1, r, g, b);
            }
            if (!visible_block(coord + glm::ivec3(0, -1, 0))) {
                verts.emplace_back(x,  y,  z,  r, g, b);
                verts.emplace_back(x1, y,  z1, r, g, b);
                verts.emplace_back(x,  y,  z1, r, g, b);
                verts.emplace_back(x,  y,  z,  r, g, b);
                verts.emplace_back(x1, y,  z,  r, g, b);
                verts.emplace_back(x1, y,  z1, r, g, b);
            }
            if (!visible_block(coord + glm::ivec3(0, 1, 0))) {
                verts.emplace_back(x1, y1, z1, r, g, b);
                verts.emplace_back(x1, y1, z,  r, g, b);
                verts.emplace_back(x,  y1, z,  r, g, b);
                verts.emplace_back(x,  y1, z1, r, g, b);
                verts.emplace_back(x1, y1, z1, r, g, b);
                verts.emplace_back(x,  y1, z,  r, g, b);
            }
            if (!visible_block(coord + glm::ivec3(0, 0, -1))) {
                verts.emplace_back(x,  y,  z,  r, g, b);
                verts.emplace_back(x,  y1, z,  r, g, b);
                verts.emplace_back(x1, y,  z,  r, g, b);
                verts.emplace_back(x1, y1, z,  r, g, b);
                verts.emplace_back(x1, y,  z,  r, g, b);
                verts.emplace_back(x,  y1, z,  r, g, b);
            }
            if (!visible_block(coord + glm::ivec3(0, 0, 1))) {
                verts.emplace_back(x,  y1, z1, r, g, b);
                verts.emplace_back(x1, y,  z1, r, g, b);
                verts.emplace_back(x1, y1, z1, r, g, b);
                verts.emplace_back(x1, y,  z1, r, g, b);
                verts.emplace_back(x,  y1, z1, r, g, b);
                verts.emplace_back(x,  y,  z1, r, g, b);
            }
        };

        for (int z = 0; z < chunk_size; ++z) {
            for (int y = 0; y < chunk_size; ++y) {
                for (int x = 0; x < chunk_size; ++x) {
                    verts_add_block(glm::ivec3(x, y, z));
                }
            }
        }
        // for (MeshVertex mv : verts) {
        //     float x = float(mv.packed_vertex & 255);
        //     float y = float((mv.packed_vertex >> 8) & 255);
        //     float z = float((mv.packed_vertex >> 16) & 255);
        //     float red   = ((mv.packed_color >> 16) & 255) * (1./255.);
        //     float green = ((mv.packed_color >> 8) & 255) * (1./255.);
        //     float blue  = (mv.packed_color & 255) * (1./255.);
        //     printf("%f %f %f   %f %f %f\n", x, y, z, red, green, blue);
        // }
        return verts;
    }

    // Return the number of bytes that will be reserved in a VBO to
    // store the given array of mesh vertices. This is not the number
    // needed, but the number that will be allocated when we have a
    // choice.
    static GLsizeiptr verts_reserved_bytes(
        const std::vector<MeshVertex>& verts)
    {
        return GLsizeiptr((64+verts.size()) * sizeof(verts[0]));
    }

    // Return the actual number of bytes needed to store the given
    // mesh vertex array.
    static GLsizeiptr verts_needed_bytes(
        const std::vector<MeshVertex>& verts)
    {
        return GLsizeiptr(verts.size() * sizeof(verts[0]));
    }

    // Fill the given MeshEntry with mesh data for the given chunk
    // group from the given world.
    static void replace(PositionedChunkGroup& pcg,
                        VoxelWorld& world,
                        MeshEntry* entry)
    {
        entry->world_id = world.id();
        entry->group_coord = group_coord(pcg);

        GLsizeiptr bytes_to_alloc = 0;
        for (int z = 0; z < edge_chunks; ++z) {
            for (int y = 0; y < edge_chunks; ++y) {
                for (int x = 0; x < edge_chunks; ++x) {
                    ChunkMesh& mesh = entry->mesh_array[z][y][x];
                    Chunk& chunk = group(pcg).chunk_array[z][y][x];
                    glm::ivec3 residue = glm::ivec3(x,y,z) * chunk_size;
                    std::vector<MeshVertex> verts =
                        make_mesh_verts(chunk, residue);
                    chunk.mesh_dirty = false;
                    auto reserved_bytes = verts_reserved_bytes(verts);

                    mesh.vbo_byte_offset = bytes_to_alloc;
                    bytes_to_alloc += reserved_bytes;
                    mesh.vbo_bytes_reserved = reserved_bytes;
                    mesh.verts = std::move(verts);
                }
            }
        }

        if (entry->vbo_name == 0 or bytes_to_alloc > entry->vbo_bytes) {
            if (entry->vbo_name == 0) glGenBuffers(1, &entry->vbo_name);
            glBindBuffer(GL_ARRAY_BUFFER, entry->vbo_name);
            glBufferData(GL_ARRAY_BUFFER, bytes_to_alloc,
                         nullptr, GL_DYNAMIC_DRAW);
            entry->vbo_bytes = bytes_to_alloc;
            entry->bytes_used = bytes_to_alloc;
            PANIC_IF_GL_ERROR;
        }

        // Parasitically depend on the MeshEntry update function.
        // This is a bit dicey -- the main violation is that MeshEntry
        // doesn't accurately describe GPU state now, as the actual
        // data has not been uploaded yet.
        update(pcg, world, entry, true);
    }

    // Update the given mesh entry with new data from dirty chunks.
    // (The always_dirty flag is used to allow replace to depend on
    // this function, if true, re-upload all chunks' data).
    //
    // Requires that the MeshEntry corresponds to the given chunk
    // group of the given world.
    static void update(PositionedChunkGroup& pcg,
                       VoxelWorld& world,
                       MeshEntry* entry,
                       bool always_dirty = false)
    {
        assert(entry->world_id == world.id());
        assert(entry->group_coord == group_coord(pcg));
        assert(entry->vbo_name != 0);

        // First step is to recompute the mesh vertices of any chunk
        // that needs it. This is also when we find out if a
        // reallocation needs to be done due to some mesh growing too
        // much.
        //
        // Don't reset the chunk dirty flag yet since we'll need it to
        // know if re-upload has to happen.
        bool requires_reallocation = false;
        GLsizeiptr bytes_to_alloc = 0; // Only used in case of realloc.
        for (int z = 0; z < edge_chunks; ++z) {
            for (int y = 0; y < edge_chunks; ++y) {
                for (int x = 0; x < edge_chunks; ++x) {
                    ChunkMesh& mesh = entry->mesh_array[z][y][x];
                    std::vector<MeshVertex>& verts = mesh.verts;
                    Chunk& chunk = group(pcg).chunk_array[z][y][x];
                    glm::ivec3 residue = glm::ivec3(x,y,z) * chunk_size;
                    bool dirty = always_dirty | chunk.mesh_dirty;

                    if (dirty)  {
                        verts = make_mesh_verts(chunk, residue);
                    }
                    auto verts_bytes = verts_needed_bytes(verts);
                    requires_reallocation |=
                        (verts_bytes > mesh.vbo_bytes_reserved);
                    bytes_to_alloc += verts_reserved_bytes(verts);
                }
            }
        }

        // Now actually upload the data. This depends on whether
        // reallocation is needed or not. We also clear the mesh dirty
        // flags at this time, now that we're actually uploading the
        // mesh data.
        if (requires_reallocation) {
            // If realloc is required, resize the VBO and re-upload
            // every chunk's mesh (I may find a more efficient way one
            // day).
            glBindBuffer(GL_ARRAY_BUFFER, entry->vbo_name);
            glBufferData(GL_ARRAY_BUFFER, bytes_to_alloc,
                         nullptr, GL_DYNAMIC_DRAW);
            entry->vbo_bytes = bytes_to_alloc;
            entry->bytes_used = bytes_to_alloc;
            PANIC_IF_GL_ERROR;

            // fprintf(stderr, "Reallocating VBO for chunk group (%i %i %i) "
            //                 "to %i bytes.\n",
            //                 int(entry->group_coord.x),
            //                 int(entry->group_coord.y),
            //                 int(entry->group_coord.z),
            //                 int(bytes_to_alloc));
            GLsizeiptr offset = 0;

            for (int z = 0; z < edge_chunks; ++z) {
                for (int y = 0; y < edge_chunks; ++y) {
                    for (int x = 0; x < edge_chunks; ++x) {
                        ChunkMesh& mesh = entry->mesh_array[z][y][x];
                        std::vector<MeshVertex>& verts = mesh.verts;
                        Chunk& chunk = group(pcg).chunk_array[z][y][x];
                        chunk.mesh_dirty = false;

                        auto bytes_reserved = verts_reserved_bytes(verts);
                        auto bytes_needed = verts_needed_bytes(verts);
                        mesh.vbo_byte_offset = offset;
                        mesh.vbo_bytes_reserved = bytes_reserved;

                        glBufferSubData(GL_ARRAY_BUFFER, offset,
                                        bytes_needed, verts.data());
                        offset += bytes_reserved;
                    }
                }
            }
            PANIC_IF_GL_ERROR;
        }
        // No reallocation case. Just re-upload the dirty chunks.
        else {
            glBindBuffer(GL_ARRAY_BUFFER, entry->vbo_name);
            PANIC_IF_GL_ERROR;
            for (int z = 0; z < edge_chunks; ++z) {
                for (int y = 0; y < edge_chunks; ++y) {
                    for (int x = 0; x < edge_chunks; ++x) {
                        ChunkMesh& mesh = entry->mesh_array[z][y][x];
                        std::vector<MeshVertex>& verts = mesh.verts;
                        Chunk& chunk = group(pcg).chunk_array[z][y][x];
                        bool dirty = always_dirty | chunk.mesh_dirty;

                        if (!dirty) continue;

                        chunk.mesh_dirty = false;
                        auto bytes_needed = verts_needed_bytes(verts);
                        assert(mesh.vbo_byte_offset + bytes_needed <= entry->vbo_bytes);
                        assert(bytes_needed <= mesh.vbo_bytes_reserved);
                        glBufferSubData(GL_ARRAY_BUFFER,
                                        mesh.vbo_byte_offset,
                                        bytes_needed,
                                        verts.data());
                    }
                }
            }
            PANIC_IF_GL_ERROR;
        }
    }

    // Render, to the current framebuffer, chunks near the camera
    // using the conventional mesh-based algorithm. TODO: honor the
    // "near" part.
    void render_world_mesh_step() noexcept
    {
        MeshStore& store = camera.get_mesh_store();
        glm::mat4 residue_vp_matrix = camera.get_residue_vp();
        glm::vec3 eye_residue;
        glm::ivec3 eye_group;
        camera.get_eye(&eye_group, &eye_residue);

        static GLuint vao = 0;
        static GLuint program_id;
        static GLint mvp_matrix_idx;

        if (vao == 0) {
            program_id = make_program(mesh_vs_source, mesh_fs_source);
            mvp_matrix_idx = glGetUniformLocation(program_id, "mvp_matrix");
            assert(mvp_matrix_idx >= 0);
            glGenVertexArrays(1, &vao);
            glBindVertexArray(vao);
            PANIC_IF_GL_ERROR;
        }

        glBindVertexArray(vao);
        glUseProgram(program_id);
        PANIC_IF_GL_ERROR;

        // TODO: Actually cull far-away chunks.
        auto draw_chunk = [&] (Chunk& chunk, ChunkMesh& mesh)
        {
            if (chunk.total_visible == 0) return;
            auto vertex_offset = mesh.vbo_byte_offset / sizeof (mesh.verts[0]);
            auto vertex_count = mesh.verts.size();
            glDrawArrays(GL_TRIANGLES, vertex_offset, vertex_count);
        };

        auto draw_group = [&] (PositionedChunkGroup& pcg)
        {
            if (group(pcg).total_visible == 0) return;
            MeshEntry* entry = store.update<Renderer>(pcg, world);

            assert(entry->vbo_name != 0);
            glBindBuffer(GL_ARRAY_BUFFER, entry->vbo_name);

            // TODO: Maybe make one VAO per MeshEntry?
            glVertexAttribIPointer(
                packed_vertex_idx,
                1,
                GL_UNSIGNED_INT,
                sizeof(MeshVertex),
                (void*) offsetof(MeshVertex, packed_vertex));
            glEnableVertexAttribArray(packed_vertex_idx);

            glVertexAttribIPointer(
                packed_color_idx,
                1,
                GL_UNSIGNED_INT,
                sizeof(MeshVertex),
                (void*) offsetof(MeshVertex, packed_color));
            glEnableVertexAttribArray(packed_color_idx);
            PANIC_IF_GL_ERROR;

            // The view matrix only takes into account the eye's
            // residue coordinate, so the model position of the group
            // actually needs to be shifted by the eye's group coord.
            glm::vec3 model_offset = glm::vec3(group_coord(pcg) - eye_group)
                                   * float(group_size);
            glm::mat4 m = glm::translate(glm::mat4(1.0f), model_offset);
            glm::mat4 mvp = residue_vp_matrix * m;
            glUniformMatrix4fv(mvp_matrix_idx, 1, 0, &mvp[0][0]);

            for (int z = 0; z < edge_chunks; ++z) {
                for (int y = 0; y < edge_chunks; ++y) {
                    for (int x = 0; x < edge_chunks; ++x) {
                        Chunk& chunk = group(pcg).chunk_array[z][y][x];
                        draw_chunk(chunk, entry->mesh_array[z][y][x]);
                        PANIC_IF_GL_ERROR;
                    }
                }
            }
        };

        // TODO: Again, need to cull chunks.
        for (PositionedChunkGroup& pcg : world.group_map) {
            draw_group(pcg);
        }
    }

    // Now time to write the AABB-raycasting renderer. Here goes...

    // Fill the given RaycastEntry with the 3D texture + AABB
    // for the given chunk group from the given world.
    static void replace(PositionedChunkGroup& pcg,
                        VoxelWorld& world,
                        RaycastEntry* entry)
    {
        entry->world_id = world.id();
        entry->group_coord = group_coord(pcg);

        if (entry->vbo_name == 0) {
            glCreateBuffers(1, &entry->vbo_name);
            auto flags = GL_DYNAMIC_STORAGE_BIT;
            auto sz = sizeof(RaycastEntry::aabb_array);
            glNamedBufferStorage(entry->vbo_name, sz, nullptr, flags);
            PANIC_IF_GL_ERROR;
        }

        if (entry->texture_name == 0) {
            glGenTextures(1, &entry->texture_name);
            GLint texture_name = entry->texture_name;

            glBindTexture(GL_TEXTURE_3D, texture_name);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
            PANIC_IF_GL_ERROR;
            glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA,
                         group_size, group_size, group_size, 0,
                         GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
            PANIC_IF_GL_ERROR;
            static_assert(sizeof(VoxelTexel) == 4,
                          "Need to change texture format to match VoxelTexel");
        }

        // Again, I parasitically depend on the update RaycastEntry
        // function using an overriding always_dirty flag.
        update(pcg, world, entry, true);
    }

    // Update the given raycast entry with new data from dirty chunks
    // (or all chunks, if always_dirty is true).
    //
    // Requires that the RaycastEntry corresponds to the given chunk
    // group of the given world.
    static void update(PositionedChunkGroup& pcg,
                       VoxelWorld& world,
                       RaycastEntry* entry,
                       bool always_dirty = false)
    {
        assert(entry->world_id == world.id());
        assert(entry->group_coord == group_coord(pcg));
        assert(entry->texture_name != 0);
        assert(entry->vbo_name != 0);

        // This is actually a hell of a lot easier than the MeshEntry.
        // I just need to visit the chunks in one pass and upload
        // dirty sub-textures and AABBs (but defer upload of AABB to
        // later to reduce calls; the vbo_dirty flag is useful here).
        bool vbo_dirty = false;
        auto texture_name = entry->texture_name;
        for (int z = 0; z < edge_chunks; ++z) {
            for (int y = 0; y < edge_chunks; ++y) {
                for (int x = 0; x < edge_chunks; ++x) {
                    Chunk& chunk = group(pcg).chunk_array[z][y][x];
                    glm::ivec3 residue = glm::ivec3(x,y,z) * chunk_size;
                    bool tex_dirty = always_dirty | chunk.texture_dirty;
                    bool aabb_dirty = always_dirty | chunk.aabb_dirty;
                    vbo_dirty |= aabb_dirty;

                    if (aabb_dirty) {
                        glm::ivec3 aabb_low, aabb_high;
                        chunk.get_aabb(&aabb_low, &aabb_high);
                        entry->aabb_array[z][y][x] =
                            PackedAABB(aabb_low + residue, aabb_high + residue);
                        // int packed_aabb_low = entry->aabb_array[z][y][x].packed_low;
                        // int packed_aabb_high = entry->aabb_array[z][y][x].packed_high;
                        // int low_x = (packed_aabb_low >> 16) & 255;
                        // int low_y = (packed_aabb_low >> 8) & 255;
                        // int low_z = packed_aabb_low & 255;
                        // int high_x = (packed_aabb_high >> 16) & 255;
                        // int high_y = (packed_aabb_high >> 8) & 255;
                        // int high_z = packed_aabb_high & 255;
                        // printf("(%i %i %i) (%i %i %i)\n",
                        //     low_x, low_y, low_z, high_x, high_y, high_z);
                    }

                    if (tex_dirty) {
                        chunk.texture_dirty = false;
                        glTextureSubImage3D(
                            texture_name, 0,
                            residue.x, residue.y, residue.z,
                            chunk_size, chunk_size, chunk_size,
                            GL_RGBA, GL_UNSIGNED_BYTE,
                            chunk.voxel_texture);
                    }
                    static_assert(sizeof(VoxelTexel) == 4,
                          "Need to change texture format to match VoxelTexel");
                }
            }
        }
        PANIC_IF_GL_ERROR;

        // Now re-upload AABBs if any changed.
        glNamedBufferSubData(entry->vbo_name, 0,
                             sizeof entry->aabb_array, entry->aabb_array);
        static_assert(sizeof entry->aabb_array >= 64,
            "Suspiciously small aabb_array, did it get replaced by a ptr?");
        PANIC_IF_GL_ERROR;
    }

    // Render, to the current framebuffer, chunks around the camera
    // using the AABB-raycast algorithm. Chunks that are near
    // enough to have been drawn using the mesh algorithm will
    // not be re-drawn.
    //
    // TODO: This last sentence is not accurate at the moment.
    void render_world_raycast_step() noexcept
    {
        RaycastStore& store = camera.get_raycast_store();
        glm::mat4 residue_vp_matrix = camera.get_residue_vp();
        glm::vec3 eye_residue;
        glm::ivec3 eye_group;
        camera.get_eye(&eye_group, &eye_residue);

        static GLuint vao = 0;
        static GLuint program_id;
        static GLuint vertex_buffer_id;
        static GLuint element_buffer_id;
        static GLint mvp_matrix_id;
        // Position of camera eye relative to the origin of the group
        // (i.e. group_size times the group coordinate).
        static GLint eye_relative_group_origin_id;
        static GLint chunk_blocks_id;
        static GLint chunk_debug_id;

        if (vao == 0) {
            program_id = make_program(raycast_vs_source, raycast_fs_source);
            mvp_matrix_id = glGetUniformLocation(program_id, "mvp_matrix");
            assert(mvp_matrix_id >= 0);
            eye_relative_group_origin_id = glGetUniformLocation(program_id,
                "eye_relative_group_origin");
            assert(eye_relative_group_origin_id >= 0);
            chunk_blocks_id = glGetUniformLocation(program_id, "chunk_blocks");
            assert(chunk_blocks_id >= 0);
            chunk_debug_id = glGetUniformLocation(program_id, "chunk_debug");
            assert(chunk_debug_id >= 0);

            glGenBuffers(1, &vertex_buffer_id);
            glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_id);
            glBufferData(
                GL_ARRAY_BUFFER, sizeof unit_box_vertices,
                unit_box_vertices, GL_STATIC_DRAW);

            glGenBuffers(1, &element_buffer_id);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer_id);
            glBufferData(
                GL_ELEMENT_ARRAY_BUFFER, sizeof unit_box_elements,
                unit_box_elements, GL_STATIC_DRAW);

            glGenVertexArrays(1, &vao);
            glBindVertexArray(vao);
        }

        glBindVertexArray(vao);
        glUseProgram(program_id);
        glUniform1i(chunk_debug_id, chunk_debug);
        glActiveTexture(GL_TEXTURE0);
        glUniform1i(chunk_blocks_id, 0);
        PANIC_IF_GL_ERROR;

        // My plan is to use instanced rendering to draw the chunk AABBs
        // of this chunk group. The base box is a 1x1x1 unit box, which
        // is stretched and repositioned in the vertex shader to the
        // true AABB.
        auto draw_group = [&] (PositionedChunkGroup& pcg)
        {
            // if (group(pcg).total_visible == 0) return;
            RaycastEntry* entry = store.update<Renderer>(pcg, world);
            assert(entry->texture_name != 0);
            assert(entry->vbo_name != 0);

            glBindTexture(GL_TEXTURE_3D, entry->texture_name);

            // The view matrix only takes into account the eye's
            // residue coordinate, so the model position of the group
            // actually needs to be shifted by the eye's group coord.
            glm::vec3 model_offset = glm::vec3(group_coord(pcg) - eye_group)
                                   * float(group_size);
            glm::mat4 m = glm::translate(glm::mat4(1.0f), model_offset);
            glm::mat4 mvp = residue_vp_matrix * m;
            glUniformMatrix4fv(mvp_matrix_id, 1, 0, &mvp[0][0]);
            PANIC_IF_GL_ERROR;

            // Similarly, the eye residue needs to be shifted by the
            // group's position.
            glm::vec3 eye_relative_group_origin = eye_residue - model_offset;
            glUniform3fv(eye_relative_group_origin_id, 1,
                &eye_relative_group_origin[0]);

            // Get the vertex attribs going, I should probably make
            // a VAO per chunk group in the RaycastEntry.

            // Verts of the unit box.
            glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_id);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer_id);
            glVertexAttribPointer(
                unit_box_vertex_idx,
                4,
                GL_FLOAT,
                false,
                sizeof(float) * 4,
                (void*)0);
            glEnableVertexAttribArray(unit_box_vertex_idx);

            // AABB residue coords (packed as integers); these are per-
            // chunk and thus instanced.
            glBindBuffer(GL_ARRAY_BUFFER, entry->vbo_name);
            glVertexAttribIPointer(
                packed_aabb_low_idx,
                1,
                GL_UNSIGNED_INT,
                sizeof(PackedAABB),
                (void*) offsetof(PackedAABB, packed_low));
            glEnableVertexAttribArray(packed_aabb_low_idx);
            glVertexAttribDivisor(packed_aabb_low_idx, 1);

            glVertexAttribIPointer(
                packed_aabb_high_idx,
                1,
                GL_UNSIGNED_INT,
                sizeof(PackedAABB),
                (void*) offsetof(PackedAABB, packed_high));
            glEnableVertexAttribArray(packed_aabb_high_idx);
            glVertexAttribDivisor(packed_aabb_high_idx, 1);
            PANIC_IF_GL_ERROR;

            // Draw all edge_chunks^3 chunks in the chunk group.
            auto instances = edge_chunks * edge_chunks * edge_chunks;
            glDrawElementsInstanced(
                GL_TRIANGLES, 36, GL_UNSIGNED_SHORT, nullptr, instances);
            PANIC_IF_GL_ERROR;
        };

        // TODO: Again, need to cull chunks.
        for (PositionedChunkGroup& pcg : world.group_map) {
            draw_group(pcg);
        }
    }
};

void render_world_mesh_step(VoxelWorld& world, Camera& camera)
{
    Renderer(world, camera).render_world_mesh_step();
}

void render_world_raycast_step(VoxelWorld& world, Camera& camera)
{
    Renderer(world, camera).render_world_raycast_step();
}

MeshStore* new_mesh_store()
{
    return new MeshStore;
}

void delete_mesh_store(MeshStore* mesh_store)
{
    delete mesh_store;
}

RaycastStore* new_raycast_store()
{
    return new RaycastStore;
}

void delete_raycast_store(RaycastStore* raycast_store)
{
    delete raycast_store;
}

void viewport(int x, int y)
{
    glViewport(0, 0, x, y);
    PANIC_IF_GL_ERROR;
};

void gl_first_time_setup()
{
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
}

void gl_clear()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
}

// Old skybox code I copied.
void load_cubemap_face(GLenum face, const char* filename)
{
    std::string full_filename = expand_filename(filename);
    SDL_Surface* surface = SDL_LoadBMP(full_filename.c_str());
    if (surface == nullptr) {
        panic(SDL_GetError(), full_filename.c_str());
    }
    if (surface->w != 1024 || surface->h != 1024) {
        panic("Expected 1024x1024 texture", full_filename.c_str());
    }
    if (surface->format->format != SDL_PIXELFORMAT_BGR24) {
        fprintf(stderr, "%i\n", (int)surface->format->format);
        panic("Expected 24-bit BGR bitmap", full_filename.c_str());
    }

    glTexImage2D(face, 0, GL_RGB, 1024, 1024, 0,
                  GL_BGR, GL_UNSIGNED_BYTE, surface->pixels);

    SDL_FreeSurface(surface);
}

GLuint load_cubemap()
{
    GLuint id = 0;
    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_CUBE_MAP, id);

    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_LOD, 0);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAX_LOD, 8);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAX_LEVEL, 8);

    load_cubemap_face(GL_TEXTURE_CUBE_MAP_NEGATIVE_X, "left.bmp");
    load_cubemap_face(GL_TEXTURE_CUBE_MAP_POSITIVE_X, "right.bmp");
    load_cubemap_face(GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, "bottom.bmp");
    load_cubemap_face(GL_TEXTURE_CUBE_MAP_POSITIVE_Y, "top.bmp");
    load_cubemap_face(GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, "back.bmp");
    load_cubemap_face(GL_TEXTURE_CUBE_MAP_POSITIVE_Z, "front.bmp");

    glGenerateMipmap(GL_TEXTURE_CUBE_MAP);

    PANIC_IF_GL_ERROR;
    glBindTexture(GL_TEXTURE_CUBE_MAP, 0);

    return id;
}

static const char skybox_vs_source[] =
"#version 330\n"
"layout(location=0) in vec3 position;\n"
"out vec3 texture_coordinate;\n"
"uniform mat4 view_matrix;\n"
"uniform mat4 proj_matrix;\n"
"void main() {\n"
    "vec4 v = view_matrix * vec4(400*position, 0.0);\n"
    "gl_Position = proj_matrix * vec4(v.xyz, 1);\n"
    "texture_coordinate = position;\n"
"}\n";

static const char skybox_fs_source[] =
"#version 330\n"
"in vec3 texture_coordinate;\n"
"out vec4 color;\n"
"uniform samplerCube cubemap;\n"
"void main() {\n"
    "vec4 c = texture(cubemap, texture_coordinate);\n"
    "c.a = 1.0;\n"
    "color = c;\n"
    "gl_FragDepth = 0.99999;\n"
"}\n";

static const float skybox_vertices[24] = {
    -1, 1, 1,
    -1, -1, 1,
    1, -1, 1,
    1, 1, 1,
    -1, 1, -1,
    -1, -1, -1,
    1, -1, -1,
    1, 1, -1,
};

static const GLushort skybox_elements[36] = {
    7, 4, 5, 7, 5, 6,
    1, 0, 3, 1, 3, 2,
    5, 1, 2, 5, 2, 6,
    4, 7, 3, 4, 3, 0,
    0, 1, 5, 0, 5, 4,
    2, 3, 7, 2, 7, 6
};

void draw_skybox(
    glm::mat4 view_matrix, glm::mat4 proj_matrix)
{
    static bool cubemap_loaded = false;
    static GLuint cubemap_texture_id;
    if (!cubemap_loaded) {
        cubemap_texture_id = load_cubemap();
        cubemap_loaded = true;
    }

    static GLuint vao = 0;
    static GLuint program_id;
    static GLuint vertex_buffer_id;
    static GLuint element_buffer_id;
    static GLint view_matrix_id;
    static GLint proj_matrix_id;
    static GLint cubemap_uniform_id;

    if (vao == 0) {
        program_id = make_program(skybox_vs_source, skybox_fs_source);
        view_matrix_id = glGetUniformLocation(program_id, "view_matrix");
        proj_matrix_id = glGetUniformLocation(program_id, "proj_matrix");
        cubemap_uniform_id = glGetUniformLocation(program_id, "cubemap");

        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        glGenBuffers(1, &vertex_buffer_id);
        glGenBuffers(1, &element_buffer_id);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer_id);
        glBufferData(
            GL_ELEMENT_ARRAY_BUFFER, sizeof skybox_elements,
            skybox_elements, GL_STATIC_DRAW
        );

        glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_id);
        glBufferData(
            GL_ARRAY_BUFFER, sizeof skybox_vertices,
            skybox_vertices, GL_STATIC_DRAW
        );
        glVertexAttribPointer(
            0,
            3,
            GL_FLOAT,
            false,
            sizeof(float) * 3,
            (void*)0
        );
        glEnableVertexAttribArray(0);
        PANIC_IF_GL_ERROR;
    }

    glUseProgram(program_id);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap_texture_id);
    glUniform1i(cubemap_uniform_id, 0);

    glUniformMatrix4fv(view_matrix_id, 1, 0, &view_matrix[0][0]);
    glUniformMatrix4fv(proj_matrix_id, 1, 0, &proj_matrix[0][0]);

    glBindVertexArray(vao);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_SHORT, (void*)0);
    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_CUBE_MAP, 0);

    PANIC_IF_GL_ERROR;
}

} // end namespace
