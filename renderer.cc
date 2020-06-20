#include "myricube.hh"

#include <stdio.h>
#include <utility>
#include <vector>

#include "camera.hh"
#include "chunk.hh"
#include "renderer.hh"
#include "glad/glad.h"

namespace myricube {

static constexpr int edge_chunks = group_size / chunk_size;

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
        packed_color = x | uint32_t(y) << 8 | uint32_t(z) << 16;
        packed_vertex = blue | uint32_t(green) << 8 | uint32_t(red) << 16;
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
        fprintf(stderr, "Evicting MeshStore entry for chunk group "
                        "(%i %i %i) of world %li\n",
                        int(entry->group_coord.x),
                        int(entry->group_coord.y),
                        int(entry->group_coord.z),
                        long(entry->world_id));

        EntryFiller::replace(pcg, world, entry);
        return entry;
    }
};

class MeshStore : public BaseStore<MeshEntry, 16> { };

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

auto packed_vertex_idx = 0;
auto packed_color_idx = 1;

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

            printf("%i %i %i\n", int(coord.x), int(coord.y), int(coord.z));

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
        return verts;
    }

    // Return the number of bytes that will be reserved in a VBO to
    // store the given array of mesh vertices. This is not the number
    // needed, but the number that will be allocated when we have a
    // choice.
    static GLsizeiptr verts_reserved_bytes(
        const std::vector<MeshVertex>& verts)
    {
        return GLsizeiptr((128+verts.size()) * sizeof(verts[0]));
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
                       bool always_dirty = false) noexcept
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

    void render_world_mesh_step()
    {
        MeshStore& store = camera.get_mesh_store();
        glm::mat4 residue_vp_matrix = camera.get_residue_vp();
        glm::vec3 eye_residue;
        glm::ivec3 eye_group;
        camera.get_eye(&eye_group, &eye_residue);

        static GLuint vao = 0;
        static GLuint program_id;
        static GLint mvp_matrix_idx = 0;

        if (vao == 0) {
            program_id = make_program(mesh_vs_source, mesh_fs_source);
            mvp_matrix_idx = glGetUniformLocation(program_id, "mvp_matrix");
            glGenVertexArrays(1, &vao);
            glBindVertexArray(vao);
            PANIC_IF_GL_ERROR;
        }

        glBindVertexArray(vao);
        glUseProgram(program_id);
        PANIC_IF_GL_ERROR;

        // TODO: Actually cull far-away chunks.
        auto draw_chunk = [&] (Chunk&, ChunkMesh& mesh)
        {
            auto vertex_offset = mesh.vbo_byte_offset / sizeof (mesh.verts[0]);
            auto vertex_count = mesh.verts.size();
            if (vertex_count != 0) {
                printf("vertex_offset %i\n", int(vertex_offset));
            }
            glDrawArrays(GL_TRIANGLES, vertex_offset, vertex_count);
        };

        auto draw_group = [&] (PositionedChunkGroup& pcg)
        {
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
};

void render_world_mesh_step(VoxelWorld& world, Camera& camera)
{
    Renderer(world, camera).render_world_mesh_step();
}

MeshStore* new_mesh_store()
{
    return new MeshStore;
}

void delete_mesh_store(MeshStore* mesh_store)
{
    delete mesh_store;
}

// TODO
RaycastStore* new_raycast_store()
{
    return nullptr;
}

void delete_raycast_store(RaycastStore*) { }

void on_window_resize(int x, int y)
{
    glViewport(0, 0, x, y);
    PANIC_IF_GL_ERROR;
};

void gl_first_time_setup()
{
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glClearColor(0.3, 0.3, 0.3, 1);
}

void gl_clear()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
}

} // end namespace
