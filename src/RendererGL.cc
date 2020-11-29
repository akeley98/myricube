// OpenGL implementation of RendererLogic.
#include "myricube.hh"
#include "EnvVar.hh"

namespace {
myricube::EnvVar64 validation_enabled("myricube_validation", 0);
}

#define PANIC_IF_GL_ERROR do { \
    if (validation_enabled) { \
        if (GLenum PANIC_error = glGetError()) { \
            char PANIC_msg[160]; \
            snprintf(PANIC_msg, sizeof PANIC_msg, "line %i: code %u", __LINE__, (unsigned)PANIC_error); \
            panic("OpenGL error", PANIC_msg); \
        } \
    } \
} while (0)

#include "glad/glad.h"
#include "RendererLogic.hh"

#include <cassert>
#include "shaders.hh"

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



// aabb_array[z][y][x] is the minimal AABB containing the visible
// voxels of ChunkGroup::chunk_array[z][y][x].
//
// If you change this, be aware that sizeof is used on this to
// know the amount of bytes to allocate for the GPU's VBO.
struct AABBs
{
    PackedAABB aabb_array[edge_chunks][edge_chunks][edge_chunks];
};

// All voxels within a chunk group share a single 3D voxel array on
// the GPU (stored as a 3D texture), as well as one VBO used to store
// the array of AABB for the chunks within the group. This is the
// handle for the GPU data needed to raycast one chunk group.
struct RaycastEntry
{
    // OpenGL name of the VBO storing the AABBs array. 0 when not yet
    // allocated.
    GLuint vbo_name = 0;

    // Memory-map view of VBO.
    AABBs* mapped_aabb = nullptr;

    // OpenGL name for the 3D voxels texture. 0 when not yet allocated.
    GLuint texture_name_ = 0;

    void bind_texture()
    {
        assert(texture_name_ != 0);
        glBindTexture(GL_TEXTURE_3D, texture_name_);
    }

    RaycastEntry()
    {
        glCreateTextures(GL_TEXTURE_3D, 1, &texture_name_);

        glBindTexture(GL_TEXTURE_3D, texture_name_);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
        PANIC_IF_GL_ERROR;
        glTexStorage3D(
            GL_TEXTURE_3D, 1, GL_RGBA8, group_size, group_size, group_size);
        PANIC_IF_GL_ERROR;

        auto storage_flags = GL_MAP_WRITE_BIT
                           | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT;
        auto map_flags = GL_MAP_PERSISTENT_BIT | GL_MAP_INVALIDATE_BUFFER_BIT
              | GL_MAP_COHERENT_BIT | GL_MAP_WRITE_BIT;

        glCreateBuffers(1, &vbo_name);
        glNamedBufferStorage(vbo_name, sizeof(AABBs),
            nullptr, storage_flags);
        mapped_aabb = static_cast<AABBs*>(glMapNamedBufferRange(
            vbo_name, 0, sizeof(AABBs), map_flags));
        PANIC_IF_GL_ERROR;
    }

    ~RaycastEntry()
    {
        GLuint buffers[1] = { vbo_name };
        glDeleteBuffers(1, buffers);
        vbo_name = 0;
        mapped_aabb = nullptr;
        glDeleteTextures(1, &texture_name_);
        texture_name_ = 0;
    }

    RaycastEntry(RaycastEntry&& other) = delete;
};

// Staging buffer for the async cache of raycast chunk groups.
struct RaycastStaging
{
   // OpenGL name of the staging SSBO.
    GLuint ssbo_name = 0;

    // Persistent mapped SSBO pointer.
    ChunkGroupVoxels* mapped_ssbo = nullptr;

    // Set to non-null once the compute shader is dispatched.
    GLsync sync = nullptr;

    // RaycastEntry to be created.
    std::unique_ptr<RaycastEntry> entry;

    // Since this will be filled by worker threads (not owning a GL
    // context), we need to make sure the buffer is 100% ready to go
    // as soon as it's created.
    RaycastStaging()
    {
        re_init();
    }

    void re_init()
    {
        auto storage_flags = GL_MAP_WRITE_BIT
                           | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT;

        auto map_flags = GL_MAP_PERSISTENT_BIT | GL_MAP_INVALIDATE_BUFFER_BIT
              | GL_MAP_COHERENT_BIT | GL_MAP_WRITE_BIT;

        if (entry == nullptr) {
            entry.reset(new RaycastEntry);
        }
        if (ssbo_name == 0) {
            glCreateBuffers(1, &ssbo_name);
            glNamedBufferStorage(ssbo_name,
                sizeof(ChunkGroupVoxels), nullptr, storage_flags);
            mapped_ssbo = static_cast<ChunkGroupVoxels*>(glMapNamedBufferRange(
                ssbo_name, 0, sizeof(ChunkGroupVoxels), map_flags));
        }
        PANIC_IF_GL_ERROR;
    }

    RaycastStaging(RaycastStaging&&) = delete;

    ~RaycastStaging()
    {
        glDeleteBuffers(1, &ssbo_name);
    }
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
        p_window->gl_make_current(gladLoadGL);

        glClearColor(0, 0, 0, 1);

        glDepthFunc(GL_LESS);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glProvokingVertex(GL_FIRST_VERTEX_CONVENTION);
        glClipControl(GL_UPPER_LEFT, GL_ZERO_TO_ONE);
        glFrontFace(GL_CW);

        // Warn if validation disabled.
        if (!validation_enabled) {
            fprintf(stderr,
            #if !defined(MYRICUBE_WINDOWS)
                "\x1b[35m\x1b[1m"
            #endif
                "Note for developers: OpenGL glGetError checks disabled!\n"
            #if !defined(MYRICUBE_WINDOWS)
                "\x1b[0m"
            #endif
                "(enable with environment variable myricube_validation=1)\n"
            );
        }
    }

    void begin_frame() override
    {
        int x, y;
        p_window->get_framebuffer_size(&x, &y);
        glViewport(0, 0, x, y);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glDisable(GL_DEPTH_TEST);
        render_background();
        glEnable(GL_DEPTH_TEST);
    }


    /* MESH RENDERING IMPLEMENTATION */

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

    // Draw the list of chunk groups with the instanced-voxel
    // rendering method.
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



    /* RAYCAST RENDERER IMPLEMENTATION */

    // Big-picture of data: I'm storing voxel data for chunk groups as
    // 3D textures. These textures live in RaycastEntry inside the
    // main cache of AsyncCache, hidden in RendererLogic.
    //
    // To fill these textures, worker threads in AsyncCache fill the
    // staging buffers (memory-mapped SSBO) with voxel data loaded
    // from disk. This is separate from the main cache.
    //
    // Once the worker thread marks a staging buffer as fully loaded,
    // it can be swapped into the main cache. This is done with a GL
    // compute shader using imageStore to transfer the data from the
    // SSBO to the texture. The main thread dispatches the shader.
    // NOTE: This is really done in two steps (two calls of
    // worker_stage, itself triggered by swap_in_from_staging in
    // RendererLogic) to give the compute shader time to finish. See
    // RaycastStaging::sync at time of writing.

    GLuint raycast_program = 0;
    GLuint raycast_vao = 0;

    struct RaycastUniformIDs
    {
        GLint mvp_matrix;
        // Position of camera eye relative to the origin of the group
        // (origin == group_size times the group coordinate).
        GLint eye_relative_group_origin;
        GLint far_plane_squared;
        GLint raycast_thresh_squared;
        GLint fog_enabled;
        GLint black_fog;
        GLint chunk_debug;
        GLint chunk_group_texture;
    } raycast_uniform_ids;

    // Draw the given list of chunk groups with the raycast-AABB
    // rendering method.
    void draw_raycast_entries(
        const std::vector<std::pair<RaycastEntry*, glm::ivec3>>& entries)
    override
    {
        auto& tr = transforms;
        glm::mat4 residue_vp_matrix = tr.residue_vp_matrix;
        glm::vec3 eye_residue = tr.eye_residue;
        glm::ivec3 eye_group = tr.eye_group;
        auto far_plane = tr.far_plane;
        auto& u = raycast_uniform_ids;

        if (raycast_vao == 0) {
            PANIC_IF_GL_ERROR;
            auto program = make_program({
                "raycast.vert", "raycast.frag",
                "fog_border.frag", "read_group_voxel.frag" });
            raycast_program = program;

            u.mvp_matrix = glGetUniformLocation(program, "mvp_matrix");
            assert(u.mvp_matrix >= 0);
            u.eye_relative_group_origin = glGetUniformLocation(program,
                "eye_relative_group_origin");
            assert(u.eye_relative_group_origin >= 0);
            u.chunk_debug = glGetUniformLocation(program, "chunk_debug");
            assert(u.chunk_debug >= 0);

            u.far_plane_squared = glGetUniformLocation(program,
                "far_plane_squared");
            assert(u.far_plane_squared >= 0);
            u.raycast_thresh_squared = glGetUniformLocation(program,
                "raycast_thresh_squared");
            assert(u.raycast_thresh_squared >= 0);
            u.fog_enabled = glGetUniformLocation(program,
                "fog_enabled");
            assert(u.fog_enabled >= 0);
            u.black_fog = glGetUniformLocation(program, "black_fog");
            assert(u.black_fog >= 0);
            u.chunk_group_texture =
                glGetUniformLocation(program, "chunk_group_texture");
            assert(u.chunk_group_texture >= 0);

            glGenVertexArrays(1, &raycast_vao);
            glBindVertexArray(raycast_vao);
        }

        glBindVertexArray(raycast_vao);
        glUseProgram(raycast_program);

        // Set uniforms unchanged per-chunk-group.
        glUniform1i(u.chunk_debug, tr.chunk_debug);
        glUniform1i(u.fog_enabled, tr.use_fog);
        glUniform1i(u.black_fog, tr.use_black_fog);
        glUniform1i(u.far_plane_squared, far_plane * far_plane);
        auto raycast_thr = raycast_threshold;
        glUniform1i(u.raycast_thresh_squared, raycast_thr * raycast_thr);

        // I'll use texture 0 for the chunk group texture.
        glActiveTexture(GL_TEXTURE0);
        glUniform1i(u.chunk_group_texture, 0);

        PANIC_IF_GL_ERROR;

        // My plan is to use instanced rendering to draw the chunk AABBs
        // of this chunk group. The base box is a 1x1x1 unit box, which
        // is stretched and repositioned in the vertex shader to the
        // true AABB.
        for (auto pair : entries) {
            RaycastEntry* entry = pair.first;
            glm::ivec3 group_coord = pair.second;
            assert(entry->vbo_name != 0);

            // The view matrix only takes into account the eye's
            // residue coordinate, so the model position of the group
            // actually needs to be shifted by the eye's group coord.
            glm::vec3 model_offset = glm::vec3(group_coord - eye_group)
                                   * float(group_size);
            glm::mat4 m = glm::translate(glm::mat4(1.0f), model_offset);
            glm::mat4 mvp = residue_vp_matrix * m;
            glUniformMatrix4fv(u.mvp_matrix, 1, 0, &mvp[0][0]);
            PANIC_IF_GL_ERROR;

            // Similarly, the eye residue needs to be shifted by the
            // group's position.
            glm::vec3 eye_relative_group_origin = eye_residue - model_offset;
            glUniform3fv(u.eye_relative_group_origin, 1,
                &eye_relative_group_origin[0]);

            // Get the instanced vertex attribs going (i.e. bind the
            // AABB residue coords, packed as integers).
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

            // Bind the relevant chunk group texture.
            entry->bind_texture();

            // Draw all edge_chunks^3 chunks in the chunk group.
            // As specified in raycast.vert, need 14 triangle strip drawn.
            auto instances = edge_chunks * edge_chunks * edge_chunks;
            glDrawArraysInstanced(
                GL_TRIANGLE_STRIP, 0, 14, instances);
            PANIC_IF_GL_ERROR;
        }

        PANIC_IF_GL_ERROR;
        glBindVertexArray(0);
    }

    // For now, just fill in the staging SSBO with data from the chunk
    // group.  This is all we can do here since there's no OpenGL
    // context for worker threads.
    void worker_stage(
        RaycastStaging* stage, const BinChunkGroup* group_ptr) override
    {
        for (int zL = 0; zL < edge_chunks; ++zL) {
        for (int yL = 0; yL < edge_chunks; ++yL) {
        for (int xL = 0; xL < edge_chunks; ++xL) {
            // Load raw voxel data chunk-by-chunk onto the SSBO.
            const BinChunk& bin_chunk = group_ptr->chunk_array[zL][yL][xL];
            auto& source_chunk = bin_chunk.voxel_array;
            auto& device_chunk =
                stage->mapped_ssbo->voxel_colors[zL][yL][xL];
            auto sz = sizeof(uint32_t) * chunk_size*chunk_size*chunk_size;
            assert(sz == sizeof(source_chunk));
            assert(sz == sizeof(device_chunk));
            memcpy(device_chunk, source_chunk, sz);

            // Also load the AABB for this chunk.
            auto chunk_residue = glm::ivec3(xL, yL, zL) * chunk_size;
            PackedAABB aabb(bin_chunk, chunk_residue);
            stage->entry->mapped_aabb->aabb_array[zL][yL][xL] = aabb;
        }
        }
        }
    }

    GLuint raycast_staging_program_id = 0;

    // Swapping into the cache is done in two steps:
    //
    // 1. Dispatch the compute shader to transfer data from SSBO in
    // the staging buffer to the texture in the cache entry.
    //
    // 2. Wait for that compute shader to finish, then swap into
    // the cache.
    bool swap_in(RaycastStaging* staging,
                 std::unique_ptr<RaycastEntry>* p_uptr_entry) override
    {
        if (raycast_staging_program_id == 0) {
            raycast_staging_program_id = make_program( { "staging.comp" } );
        }

        // Step 1
        if (staging->sync == nullptr) {
            glUseProgram(raycast_staging_program_id);
            // Hook up the SSBO binding to the SSBO index that will
            // be used to deliver voxel data.
            glShaderStorageBlockBinding(
                raycast_staging_program_id,
                chunk_group_voxels_program_index,
                chunk_group_voxels_binding_index);
            // Hook up the out_image to the correct image unit.
            glUniform1i(staging_image_program_index, staging_image_unit);

            // Bind the correct texture image.
            RaycastEntry& entry = *staging->entry;

            glBindBufferBase(
                GL_SHADER_STORAGE_BUFFER,
                chunk_group_voxels_binding_index,
                staging->ssbo_name);
            glBindImageTexture(
                staging_image_unit,
                entry.texture_name_,
                0, false, 0, GL_WRITE_ONLY, GL_RGBA8UI);

            glDispatchCompute(1, 2, 2);
            staging->sync = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
            assert(staging->sync != nullptr);
            return false; // Not yet done.
        }
        // Step 2
        else {
            // Wait for compute shader to finish.
            // I *think* I need the client wait too because the client
            // writes to the staging SSBO.
            glWaitSync(staging->sync, 0, GL_TIMEOUT_IGNORED);
            glClientWaitSync(staging->sync, GL_SYNC_FLUSH_COMMANDS_BIT, 1e9);
            glDeleteSync(staging->sync);
            glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT
                          | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
            staging->sync = nullptr;

            // Ready. Swap into the cache.
            std::swap(staging->entry, *p_uptr_entry);
            staging->re_init();
            return true;
        }
        assert(0);
    }



    /* BACKGROUND RENDERING */
    GLuint background_vao = 0;
    GLuint background_program = 0;

    struct BackgroundUniformIDs
    {
        GLint inverse_vp;
        GLint eye_world_position;
        GLint black_fog;
    } background_uniform_ids;

    // Render the background. This uses a shader hard-wired to draw a
    // full-screen rectangle.
    void render_background()
    {
        auto& u = background_uniform_ids;

        if (background_vao == 0) {
            glGenVertexArrays(1, &background_vao);
            background_program = make_program(
                { "background.vert", "background.frag", "fog_border.frag" } );

            u.inverse_vp = glGetUniformLocation(background_program,
                "inverse_vp");
            assert(u.inverse_vp >= 0);
            u.eye_world_position = glGetUniformLocation(background_program,
                "eye_world_position");
            assert(u.eye_world_position >= 0);
            u.black_fog = glGetUniformLocation(background_program,
                "black_fog");
            assert(u.black_fog >= 0);
            PANIC_IF_GL_ERROR;
        }

        auto& tr = transforms;
        glm::vec3 eye_residue = tr.eye_residue;
        glm::mat4 inverse_vp = glm::inverse(tr.residue_vp_matrix);

        glBindVertexArray(background_vao);
        glUseProgram(background_program);
        glUniformMatrix4fv(u.inverse_vp, 1, 0, &inverse_vp[0][0]);
        glUniform3fv(u.eye_world_position, 1, &eye_residue[0]);
        glUniform1i(u.black_fog, tr.use_black_fog);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glBindVertexArray(0);
        PANIC_IF_GL_ERROR;
    };



    void end_frame() override
    {
        p_window->swap_buffers();
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
