// OpenGL implementation of RendererLogic.
#include "RendererLogic.hh"

#include "glad/glad.h"
#include "shaders.hh"
#include <stdio.h>

using namespace myricube;

namespace {

// Handle for GPU resources for the mesh of one chunk group.
struct MeshEntry
{
    // OpenGL name for the VBO used to store the meshes of this chunk group.
    GLuint vbo_name = 0;

    // Persistent coherent mapping of the VBO.
    MappedGroupMesh* vbo_map;

    // Data needed to draw each chunk [z][y][x].
    ChunkDrawData draw_data[edge_chunks][edge_chunks][edge_chunks];

    // Size in bytes of the vbo's data store on the GPU.
    static constexpr GLsizeiptr vbo_bytes = sizeof(MappedGroupMesh);

    // Return the offset (in number of MeshVoxelVertex's, not bytes)
    // into the VBO where the data from mesh_array[z][y][x] is copied
    // into.
    static unsigned vert_offset(unsigned x, unsigned y, unsigned z)
    {
        assert(x < edge_chunks and y < edge_chunks and z < edge_chunks);
        unsigned chunk_idx = x + y*edge_chunks + z*edge_chunks*edge_chunks;
        return chunk_idx * chunk_max_verts;
    }

    // Same as above, but return offset as count of bytes.
    static GLsizeiptr byte_offset(unsigned x, unsigned y, unsigned z)
    {
        auto vert_sz = GLsizeiptr(sizeof(MeshVoxelVertex));
        GLsizeiptr off = vert_offset(x, y, z) * vert_sz;
        assert(size_t(off) < vbo_bytes);
        return off;
    }

    MeshEntry()
    {
        auto storage_flags = GL_MAP_WRITE_BIT
                           | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT;

        auto map_flags = GL_MAP_PERSISTENT_BIT | GL_MAP_INVALIDATE_BUFFER_BIT
              | GL_MAP_COHERENT_BIT | GL_MAP_WRITE_BIT;


        glCreateBuffers(1, &vbo_name);
        glNamedBufferStorage(vbo_name, sizeof(MappedGroupMesh),
            nullptr, storage_flags);
        vbo_map = static_cast<MappedGroupMesh*>(glMapNamedBufferRange(
                vbo_name, 0, sizeof(MappedGroupMesh), map_flags));
        PANIC_IF_GL_ERROR;
    }

    // Get some RAII going.
    ~MeshEntry()
    {
        glDeleteBuffers(1, &vbo_name);
        vbo_name = 0;
    }

    MeshEntry(MeshEntry&& other) = delete;
};

void swap(MeshEntry& left, MeshEntry& right) noexcept
{
    using std::swap;
    swap(left.vbo_name, right.vbo_name);
    swap(left.vbo_map, right.vbo_map);
    swap(left.draw_data, right.draw_data);
}

using MeshStaging = MeshEntry;

struct RaycastEntry
{

};

struct RaycastStaging
{

};

} // end anonymous namespace.

namespace myricube {

struct RendererGL :
    RendererLogic<MeshEntry, MeshStaging, RaycastEntry, RaycastStaging>
{
    RendererGL(RenderThread* thread, RenderArgs args) :
        RendererLogic<MeshEntry, MeshStaging, RaycastEntry, RaycastStaging>
        (thread, args)
    {
        p_window->gl_make_current();

        glClearColor(0, 0, 0, 1);

        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glProvokingVertex(GL_FIRST_VERTEX_CONVENTION);
        glClipControl(GL_UPPER_LEFT, GL_ZERO_TO_ONE);
        glFrontFace(GL_CW);
    }

    void begin_frame() override
    {
        int x, y;
        p_window->get_framebuffer_size(&x, &y);
        glViewport(0, 0, x, y);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }



    // GL identifiers needed for mesh rendering.
    GLuint mesh_vao = 0;
    GLuint mesh_program = 0;

    struct MeshUniformIDs {
        GLint mvp_matrix;

        // Position of camera eye relative to the origin of the group
        // (origin == group_size times the group coordinate).
        GLint eye_relative_group_origin;

        GLint far_plane_squared;
        GLint fog_enabled;
        GLint black_fog;
    } mesh_uniform_ids;

    void draw_mesh_entries(
        const std::vector<std::pair<MeshEntry*, glm::ivec3>>& entries) override
    {
        auto& tr = transforms;
        auto& u = mesh_uniform_ids;
        glm::mat4 residue_vp_matrix = tr.residue_vp_matrix;
        glm::vec3 eye_residue = tr.eye_residue;
        glm::ivec3 eye_group = tr.eye_group;

        if (mesh_vao == 0) {
            mesh_program = make_program({
                "mesh.vert", "mesh.frag", "fog_border.frag" });
            u.mvp_matrix = glGetUniformLocation(mesh_program, "mvp_matrix");
            assert(u.mvp_matrix >= 0);
            u.eye_relative_group_origin = glGetUniformLocation(mesh_program,
                "eye_relative_group_origin");
            assert(u.eye_relative_group_origin >= 0);
            u.far_plane_squared = glGetUniformLocation(mesh_program,
                "far_plane_squared");
            assert(u.far_plane_squared >= 0);
            u.fog_enabled = glGetUniformLocation(mesh_program,
                "fog_enabled");
            assert(u.fog_enabled >= 0);
            u.black_fog = glGetUniformLocation(mesh_program, "black_fog");
            assert(u.black_fog >= 0);

            glGenVertexArrays(1, &mesh_vao);
            glBindVertexArray(mesh_vao);
            PANIC_IF_GL_ERROR;
        }

        glBindVertexArray(mesh_vao);

        glUseProgram(mesh_program);
        auto far_plane = tr.far_plane;
        glUniform1i(u.far_plane_squared, far_plane * far_plane);
        glUniform1i(u.fog_enabled, tr.use_fog);
        glUniform1i(u.black_fog, tr.use_black_fog);
        PANIC_IF_GL_ERROR;
        unsigned drawn_group_count = 0;
        unsigned drawn_chunk_count = 0;

        auto draw_group = [&] (const MeshEntry& entry, glm::ivec3 group_coord)
        {
            assert(entry.vbo_name != 0);
            glBindBuffer(GL_ARRAY_BUFFER, entry.vbo_name);

            // Bind instanced (per-voxel) attributes.
            glVertexAttribIPointer(
                packed_vertex_idx,
                1,
                GL_UNSIGNED_INT,
                sizeof(MeshVoxelVertex),
                (void*) offsetof(MeshVoxelVertex, packed_residue_face_bits));
            glEnableVertexAttribArray(packed_vertex_idx);
            glVertexAttribDivisor(packed_vertex_idx, 1);

            glVertexAttribIPointer(
                packed_color_idx,
                1,
                GL_UNSIGNED_INT,
                sizeof(MeshVoxelVertex),
                (void*) offsetof(MeshVoxelVertex, packed_color));
            glEnableVertexAttribArray(packed_color_idx);
            glVertexAttribDivisor(packed_color_idx, 1);
            PANIC_IF_GL_ERROR;

            // The view matrix only takes into account the eye's
            // residue coordinate, so the model position of the group
            // actually needs to be shifted by the eye's group coord.
            glm::vec3 model_offset = glm::vec3(group_coord - eye_group)
                                   * float(group_size);
            glm::mat4 m = glm::translate(glm::mat4(1.0f), model_offset);
            glm::mat4 mvp = residue_vp_matrix * m;
            glUniformMatrix4fv(u.mvp_matrix, 1, 0, &mvp[0][0]);

            // Similarly, the eye residue needs to be shifted by the
            // group's position.
            glm::vec3 eye_relative_group_origin = eye_residue - model_offset;
            glUniform3fv(u.eye_relative_group_origin, 1,
                &eye_relative_group_origin[0]);

            for (int z = 0; z < edge_chunks; ++z) {
                for (int y = 0; y < edge_chunks; ++y) {
                    for (int x = 0; x < edge_chunks; ++x) {
                        auto draw_data = entry.draw_data[z][y][x];
                        PackedAABB aabb = draw_data.aabb;
                        if (decide_chunk(group_coord, aabb) != draw_mesh) {
                            continue;
                        }

                        if (draw_data.vert_count == 0) continue;

                        // Need to choose the base instance because
                        // every chunk's data is at a different offset
                        // within the VBO for instanced data (position/color).
                        //
                        // As requested by mesh.vert (includes
                        // built-in model of a cube), need to draw as
                        // GL_TRIANGLES with 36 vertices.
                        glDrawArraysInstancedBaseInstance(
                            GL_TRIANGLES,
                            0, 36,
                            draw_data.vert_count,
                            entry.vert_offset(x, y, z));
                            // ^^^ Instance offset, depends on chunk.
                        ++drawn_chunk_count;
                    }
                }
            }
            ++drawn_group_count;
            PANIC_IF_GL_ERROR;
        };

        for (auto pair : entries) {
            draw_group(*pair.first, pair.second);
        }

        glBindVertexArray(0);
    }

    // Convert chunk group's chunks to meshes and pack into the
    // memory-mapped array of chunks.
    void worker_stage(
        MeshStaging* staging, const BinChunkGroup* group_ptr) override
    {
        for (int zL = 0; zL < edge_chunks; ++zL) {
        for (int yL = 0; yL < edge_chunks; ++yL) {
        for (int xL = 0; xL < edge_chunks; ++xL) {
            fill_chunk_mesh(&staging->vbo_map->chunks[zL][yL][xL],
                            &staging->draw_data[zL][yL][xL],
                            group_ptr->chunk_array[zL][yL][xL],
                            glm::ivec3(xL, yL, zL) * chunk_size);
        }
        }
        }
    }

    // Since MeshEntry and MeshStaging are actually the same,
    // literally swap-in the staged buffer into the cache (with the
    // evicted cache entry, which is no longer used for rendering,
    // used as the new staging buffer).
    bool swap_in(
        MeshStaging* staging,
        std::unique_ptr<MeshEntry>* p_uptr_entry) override
    {
        if (*p_uptr_entry == nullptr) p_uptr_entry->reset(new MeshEntry);

        swap(*staging, **p_uptr_entry);
        return true;
    }

    void
    draw_raycast_entries(
        const std::vector<std::pair<RaycastEntry*, glm::ivec3>>&) override
    {

    }

    void end_frame() override
    {
        p_window->swap_buffers();
    }



    void worker_stage(RaycastStaging*, const BinChunkGroup*) override
    {

    }



    bool swap_in(RaycastStaging*, std::unique_ptr<RaycastEntry>*) override
    {
        return true;
    }

    void wait_idle() override
    {
        glFinish();
    }
};

std::shared_ptr<RendererBase> RendererGL_Factory(
    RenderThread* thread,
    RenderArgs args)
{
    return std::shared_ptr<RendererBase>(new RendererGL(thread, args));
}

} // end namespace myricube
