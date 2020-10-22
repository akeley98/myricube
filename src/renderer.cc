// Implementation of the hybrid voxel renderer. I don't see how I can
// make much of this exception safe, so I wrap everything in a
// noexcept at the end. In theory this is threadsafe (i.e. multiple
// Renderers drawing different voxel worlds using multiple GL
// contexts) but I haven't tested it so it almost certainly is not.

#include "myricube.hh"

#include <algorithm>
#include <atomic>
#include <memory>
#include <stdio.h>
#include <string.h>
#include <thread>
#include <typeinfo>
#include <utility>
#include <vector>

#include "AsyncCache.hh"
#include "camera.hh"
#include "glad/glad.h"
#include "renderer.hh"
#include "shaders.hh"
#include "voxels.hh"
#include "window.hh"

namespace myricube {

// Hacky debug variables.
bool chunk_debug = false;
bool evict_stats_debug = false; // Unused for now
bool disable_zcull_sort = false;

// (roughly) minimum distance from the camera that a chunk needs
// to be to switch from mesh to raycast graphics.
// Keep as int to avoid rounding errors in distance culling.
constexpr int raycast_threshold = 170;

// Slightly higher threshold than raycast_threshold. Chunk groups
// within this distance to the camera have their meshes loaded even
// if no part is actually drawn as a mesh (this is needed now that
// MeshStore has some latency in loading meshes).
constexpr int mesh_load_threshold = 10 + raycast_threshold;

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
        packed_low = uint32_t(x) << x_shift
                   | uint32_t(y) << y_shift
                   | uint32_t(z) << z_shift;
        x = aabb_high.x;
        y = aabb_high.y;
        z = aabb_high.z;
        packed_high = uint32_t(x) << x_shift
                    | uint32_t(y) << y_shift
                    | uint32_t(z) << z_shift;
    }

    // Compute the AABB (in residue coordinates) of the given chunk.
    // The residue coordinates of the lower corner of the chunk must
    // also be given (as the AABB is relative to the chunk group's
    // origin, not the chunk itself).
    PackedAABB(const BinChunk& chunk, glm::ivec3 chunk_residue)
    {
        // {xyz}_array[n] has visible_bit set iff there's some visible
        // voxel in the chunk with n == the voxel's {xyz} coordinate.
        uint32_t x_array[chunk_size] = { 0 };
        uint32_t y_array[chunk_size] = { 0 };
        uint32_t z_array[chunk_size] = { 0 };

        uint32_t all_orrd = 0;

        for (size_t z = 0; z < chunk_size; ++z) {
            for (size_t y = 0; y < chunk_size; ++y) {
                for (size_t x = 0; x < chunk_size; ++x) {
                    uint32_t voxel = chunk.voxel_array[z][y][x];
                    all_orrd |= voxel;
                    x_array[x] |= voxel;
                    y_array[y] |= voxel;
                    z_array[z] |= voxel;
                }
            }
        }

        // Return trivial PackedAABB if chunk is empty.
        if ((all_orrd & visible_bit) == 0) {
            packed_low = 0;
            packed_high = 0;
            return;
        }

        // Otherwise, we actually have to calculate the AABB. Do it
        // in the coordinate system of the chunk first.
        glm::ivec3 chunk_aabb_low(chunk_size, chunk_size, chunk_size);
        glm::ivec3 chunk_aabb_high(0, 0, 0);

        for (int32_t i = 0; i < chunk_size; ++i) {
            if (x_array[i] & visible_bit) {
                chunk_aabb_low.x = glm::min(chunk_aabb_low.x, i);
                chunk_aabb_high.x = glm::max(chunk_aabb_low.x, i+1);
            }
            if (y_array[i] & visible_bit) {
                chunk_aabb_low.y = glm::min(chunk_aabb_low.y, i);
                chunk_aabb_high.y = glm::max(chunk_aabb_low.y, i+1);
            }
            if (z_array[i] & visible_bit) {
                chunk_aabb_low.z = glm::min(chunk_aabb_low.z, i);
                chunk_aabb_high.z = glm::max(chunk_aabb_low.z, i+1);
            }
        }

        // Reposition computed AABB by the chunk's residue coordinate.
        *this = PackedAABB(
            chunk_aabb_low + chunk_residue,
            chunk_aabb_high + chunk_residue);
    }

    uint8_t low_x() const
    {
        return uint8_t((packed_low >> x_shift) & 255);
    }
    uint8_t low_y() const
    {
        return uint8_t((packed_low >> y_shift) & 255);
    }
    uint8_t low_z() const
    {
        return uint8_t((packed_low >> z_shift) & 255);
    }
    uint8_t high_x() const
    {
        return uint8_t((packed_high >> x_shift) & 255);
    }
    uint8_t high_y() const
    {
        return uint8_t((packed_high >> y_shift) & 255);
    }
    uint8_t high_z() const
    {
        return uint8_t((packed_high >> z_shift) & 255);
    }
};

// My new plan for the GPU voxel mesh: each visible voxel will be
// represented by ONE vertex in a VBO; instanced rendering will
// transform this vertex into an actual cube. This vertex encodes, in
// a packed bitfield format,
//
// The residue coordinates of the voxel
//
// The color of the voxel (8-bits each for R, G, B)
//
// Which faces of the voxel are visible (+/- x, y, z faces).
// (This depends on whether the neighboring voxels are filled or not)
struct MeshVoxelVertex
{
    uint32_t packed_residue_face_bits = 0xFFFFFFFF;
    uint32_t packed_color = 0xFFFFFFFF;

    MeshVoxelVertex() = default;

    // Given a _visible_ voxel and its residue coordinate, return the
    // VBO vertex carrying this information.
    //
    // This constructor defines the packed layout.
    MeshVoxelVertex(uint32_t packed_color_arg, uint8_t x, uint8_t y, uint8_t z)
    {
        static_assert(group_size <= 255,
                      "group too big for 8-bit unsigned coordinates.");
        packed_residue_face_bits = uint32_t(x) << x_shift
                                 | uint32_t(y) << y_shift
                                 | uint32_t(z) << z_shift;
        assert(packed_color_arg & visible_bit);
        packed_color = packed_color_arg;
    }
};

// The maximum number of MeshVoxelVerts needed for one chunk.  For now
// it's just the total number of voxels per chunk. I'm sure there's a
// lower possible bound but for now I'll be conservative (even though
// I'm desperate for GPU memory).
constexpr size_t chunk_max_verts = chunk_size * chunk_size * chunk_size;

// Bytes on GPU for storing the mesh (list of visible voxels) of one
// chunk.
struct MappedChunkMesh
{
    MeshVoxelVertex verts[chunk_max_verts];
};

// An entire chunk group of MappedChunkMesh, [z][y][x] order as typical.
struct MappedGroupMesh
{
    MappedChunkMesh chunks[edge_chunks][edge_chunks][edge_chunks];
};

// Extra data needed to interpret MappedChunkMesh correctly for drawing.
struct ChunkDrawData
{
    size_t vert_count = 0; // Vertex (visible voxel) count i.e. # instances
    PackedAABB aabb;       // AABB, for decide_chunk's benefit.
};

// Function for filling the above structures given a chunk. Requires
// in addition the residue coordinate of the chunk's lower-left corner
// (needed for positioning the AABB and voxel positions in the
// containing chunk group's coordinate system).
void fill_chunk_mesh(
    MappedChunkMesh* mesh_ptr,
    ChunkDrawData* draw_data_ptr,
    const BinChunk& chunk,
    glm::ivec3 chunk_residue)
{
    draw_data_ptr->aabb = PackedAABB(chunk, chunk_residue);

    // Look up whether the voxel at the given coordinate
    // (relative to the lower-left of this chunk) is visible.
    // Act as if voxels outside the chunk are always invisible.
    auto visible_block = [&chunk] (glm::ivec3 coord) -> bool
    {
        if (coord.x < 0 or coord.x >= chunk_size
         or coord.y < 0 or coord.y >= chunk_size
         or coord.z < 0 or coord.z >= chunk_size) return false;
        // Note: the masking in Chunk::operator () won't mess up
        // this coord. (I'm kind of violating my own comment in
        // Chunk::operator() because coord won't actually be in
        // the chunk unless that chunk is at (0,0,0)).
        return chunk(coord) & visible_bit;
    };

    auto visit_voxel = [
        mesh_ptr, draw_data_ptr, &chunk, visible_block, chunk_residue]
    (glm::ivec3 coord)
    {
        auto v = chunk(coord);
        if (0 == (v & visible_bit)) return;

        uint8_t x = uint8_t(coord.x) + chunk_residue.x;
        uint8_t y = uint8_t(coord.y) + chunk_residue.y;
        uint8_t z = uint8_t(coord.z) + chunk_residue.z;

        MeshVoxelVertex vert(v, x, y, z);

        // Check which of the six faces are visible.
        if (!visible_block(coord + glm::ivec3(-1, 0, 0))) {
            vert.packed_residue_face_bits |= neg_x_face_bit;
        }
        if (!visible_block(coord + glm::ivec3(1, 0, 0))) {
            vert.packed_residue_face_bits |= pos_x_face_bit;
        }
        if (!visible_block(coord + glm::ivec3(0, -1, 0))) {
            vert.packed_residue_face_bits |= neg_y_face_bit;
        }
        if (!visible_block(coord + glm::ivec3(0, 1, 0))) {
            vert.packed_residue_face_bits |= pos_y_face_bit;
        }
        if (!visible_block(coord + glm::ivec3(0, 0, -1))) {
            vert.packed_residue_face_bits |= neg_z_face_bit;
        }
        if (!visible_block(coord + glm::ivec3(0, 0, 1))) {
            vert.packed_residue_face_bits |= pos_z_face_bit;
        }

        // Add this voxel only if it's visible.
        if ((vert.packed_residue_face_bits & all_face_bits) != 0) {
            assert(draw_data_ptr->vert_count < chunk_max_verts);
            mesh_ptr->verts[draw_data_ptr->vert_count++] = vert;
        }
    };

    draw_data_ptr->vert_count = 0;

    for (int z = 0; z < chunk_size; ++z) {
        for (int y = 0; y < chunk_size; ++y) {
            for (int x = 0; x < chunk_size; ++x) {
                visit_voxel(glm::ivec3(x, y, z));
            }
        }
    }
}

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

// Renderer class. Stores state that is used by a voxel rendering
// thread running for the lifetime of the Renderer.
class Renderer
{
    std::shared_ptr<Window> window_ptr;
    std::shared_ptr<SyncCamera> sync_camera_ptr;
    WorldHandle world_handle;
    ViewWorldCache world_cache;
    const BinChunkGroup* current_chunk_group = nullptr;

    // Updated from SyncCamera per frame.
    CameraTransforms tr;

    // Device storage for mesh and raycast rendering respectively.
    // These need to be unique_ptr both because of incomplete type
    // errors and also since only the renderer thread (which is NOT
    // the one calling the Renderer constructor) can intialize these
    // properly.
    class MeshStore;
    class RaycastStore;
    std::unique_ptr<MeshStore> mesh_store;
    std::unique_ptr<RaycastStore> raycast_store;

    // Used for communicating with the render thread.
    std::atomic<bool> thread_exit_flag { false };
    std::atomic<double> fps { 0 };
    std::atomic<double> frame_time { 0 };

    std::thread thread;

  public:
    Renderer(
        std::shared_ptr<Window> window_,
        const WorldHandle& world_,
        std::shared_ptr<SyncCamera> sync_camera_) :

        window_ptr(std::move(window_)),
        sync_camera_ptr(std::move(sync_camera_)),
        world_handle(world_),
        world_cache(world_),
        thread()
    {
        thread = std::thread(render_loop, this);
    }

    ~Renderer()
    {
        thread_exit_flag.store(true);
        thread.join();
    }

    Renderer(Renderer&&) = delete;

    friend double get_fps(const Renderer& renderer);
    friend double get_frame_time(const Renderer& renderer);

  private:
    // Sets the current_chunk_group ptr to point to the named chunk group.
    // (nullptr if it doesn't exist). This seems oddly stateful for my
    // tastes but I did it anyway to avoid causing a dangling pointer.
    // (Subsequent world_cache.get_entry may overwrite ptr).
    //
    // I may re-think the design of this program at some point. It used to
    // be quite clean but the quality has eroded with time as I stapled-on
    // more optimizations. But for now I'm just having fun.
    void set_current_chunk_group(glm::ivec3 group_coord)
    {
        auto& entry = world_cache.get_entry(group_coord);
        current_chunk_group = entry.chunk_group_ptr.get();
    }

    // Culling strategy: There are two phases to culling chunks: group
    // culling and distance culling. In group culling, entire groups of
    // chunks are culled if they are outside the view frustum, or if they
    // contain no visible voxels. Then of the chunks in non-culled groups,
    // their distance to the eye is calculated, and the chunk is culled,
    // drawn with raycasting, or drawn as a mesh, depending on whether
    // it is out-of-range of, far, or near the camera.

    // Should this group be culled? If not, return through
    // *squared_dist the squared distance between the eye and the
    // nearest point in the group.
    bool cull_group(glm::ivec3 group_coord, float* squared_dist=nullptr)
    {
        glm::mat4 vp = tr.residue_vp_matrix;
        glm::ivec3 eye_group = tr.eye_group;
        glm::vec3 eye_residue = tr.eye_residue;

        // Never cull the group the eye is in (this avoids
        // pathological cases e.g. 7 corners behind eye and 1 in front
        // and out of the frustum).
        if (group_coord == eye_group) {
            if (squared_dist) *squared_dist = 0.0f;
            return false;
        }

        eye_residue = glm::floor(eye_residue); // to match decide_chunk.

        // Position of this chunk group relative to the group that the
        // eye is in, in voxel units.
        glm::vec3 low_corner(group_size * (group_coord - eye_group));

        const glm::vec3 x_edge = glm::vec3(group_size, 0, 0);
        const glm::vec3 y_edge = glm::vec3(0, group_size, 0);
        const glm::vec3 z_edge = glm::vec3(0, 0, group_size);

        // Compute the clip space (?) coordinates of the 8 corners of
        // this chunk group, and find the min/max of the xyz coordinates.
        // Take the abs of w before dividing x and y (but not z) so that
        // the x/y bounding planes don't swap sides.
        bool first_time = true;
        glm::vec3 low(-1), high(-1);
        auto minmax_corner = [vp, &low, &high, &first_time] (glm::vec3 v)
        {
            auto vp_v = vp * glm::vec4(v, 1);
            float abs_w = glm::abs(vp_v.w);
            glm::vec3 clip_coord(vp_v.x/abs_w, vp_v.y/abs_w, vp_v.z/vp_v.w);

            if (first_time) {
                low = clip_coord;
                high = clip_coord;
                first_time = false;
            }
            else {
                low = glm::min(low, clip_coord);
                high = glm::max(high, clip_coord);
            }
        };

        minmax_corner(low_corner);
        minmax_corner(low_corner + x_edge);
        minmax_corner(low_corner + y_edge);
        minmax_corner(low_corner + z_edge);
        minmax_corner(low_corner + x_edge + y_edge);
        minmax_corner(low_corner + x_edge + z_edge);
        minmax_corner(low_corner + y_edge + z_edge);
        glm::vec3 high_corner = low_corner + x_edge + y_edge + z_edge;
        minmax_corner(high_corner);

        // Cull the chunk if it is entirely out-of-bounds in one
        // direction on the x/y/z axis (one-direction == we won't clip
        // the group if the corners are all out-of-bounds but some
        // visible portion of the group "straddles" the frustum).
        if (low.x < -1 && high.x < -1) return true;
        if (low.y < -1 && high.y < -1) return true;
        if (low.z < 0 && high.z < 0) return true;
        if (low.x > 1 && high.x > 1) return true;
        if (low.y > 1 && high.y > 1) return true;
        if (low.z > 1 && high.z > 1) return true;

        if (squared_dist) {
            glm::vec3 nearest_point = glm::clamp(
                eye_residue, low_corner, high_corner);
            glm::vec3 disp = nearest_point - eye_residue;
            *squared_dist = glm::dot(disp, disp);
        }
        return false;
    }

    static constexpr int cull = 1, draw_mesh = 2, draw_raycast = 3;
    // Distance culling step. Given the group coordinate and AABB of a
    // chunk, decide which of the above things we should do: cull it,
    // draw it with the mesh algorithm (near the eye), or draw it with
    // the raycast algorithm (far away).
    //
    // This function may be duplicated on the GPU, hence the
    // floor(eye) -- this prevents subtle disagreements due to
    // different FP representations.
    int decide_chunk(glm::ivec3 group_coord,
                     PackedAABB aabb)
    {
        glm::ivec3 aabb_low(aabb.low_x(), aabb.low_y(), aabb.low_z());
        glm::ivec3 aabb_high(aabb.high_x(), aabb.high_y(), aabb.high_z());

        glm::ivec3 eye_group = tr.eye_group;
        glm::vec3 eye_residue = tr.eye_residue;
        auto far_plane = tr.far_plane;
        auto raycast_thresh = raycast_threshold;

        // Compute the displacement between the eye and the nearest point
        // of the chunk's AABB (to the eye).
        glm::ivec3 floor_eye = glm::ivec3(glm::floor(eye_residue))
                             + group_size * (eye_group - group_coord);
        glm::ivec3 aabb_nearest = glm::clamp(floor_eye, aabb_low, aabb_high);
        auto disp = glm::vec3(aabb_nearest - floor_eye);

        float squared_dist = glm::dot(disp, disp);
        if (squared_dist > far_plane * far_plane) return cull;
        if (squared_dist < raycast_thresh * raycast_thresh) return draw_mesh;
        return draw_raycast;
    }

    // Storage class for chunk group meshes.
    class MeshStore : public AsyncCache<MeshEntry, MeshEntry>
    {
        // Back pointer to the Renderer.
        Renderer* owner;

        static constexpr AsyncCacheArgs args = {
            3,          // modulus
            8,          // associativity
            32,         // staging buffers
            1,          // worker threads
            10,         // condvar timeout in milliseconds
        };

      public:
        MeshStore(Renderer* owner_) :
            AsyncCache(args), owner(owner_)
        {

        }

        // Read from disk/memory data for the given chunk group and
        // fill in the mesh.
        bool stage(MeshEntry* staging, glm::ivec3 group_coord) override
        {
            thread_local std::unique_ptr<ViewWorldCache> world_cache_ptr;

            if (world_cache_ptr == nullptr) {
                world_cache_ptr.reset(new ViewWorldCache(owner->world_handle));
            }

            auto& entry = world_cache_ptr->get_entry(group_coord);
            const BinChunkGroup* group_ptr = entry.chunk_group_ptr.get();
            if (group_ptr == nullptr) return false;

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
            return true;
        }

        // Swap from staging buffer into the cache. Since they're both
        // the same type, I just use the swap function. NOTE: I need
        // to figure out a way to prevent still-in-use entries from
        // being swapped out.
        bool swap_in(MeshEntry* staging,
                     std::unique_ptr<MeshEntry>* p_uptr_entry) override
        {
            if (*p_uptr_entry == nullptr) p_uptr_entry->reset(new MeshEntry);

            swap(**p_uptr_entry, *staging);
            return true;
        }
    };

    // GL objects needed for mesh rendering.
    GLuint mesh_vao = 0;
    GLuint mesh_program = 0;

    // Render, to the current framebuffer, chunks near the camera
    // using the conventional mesh-based algorithm.
    //
    // Each MeshEntry contains a list of visible voxels: their residue
    // coordinates, colors, and which of their 6 faces are visible.
    // My plan is to use instanced rendering to draw a bunch of unit
    // cubes in the correct locations, using degenerate triangles to
    // hide the hidden faces.
    //
    // Here we just collect a list of MeshEntries near the camera to
    // draw.  We shunt off the real work to draw_mesh_entries.
    void render_world_mesh_step() noexcept
    {
        std::vector<std::pair<MeshEntry*, glm::ivec3>> entries;
        MeshStore& store = *mesh_store;
        store.begin_frame();

        float squared_thresh = mesh_load_threshold * mesh_load_threshold;

        glm::ivec3 group_coord_low, group_coord_high;
        glm::dvec3 disp(mesh_load_threshold);
        split_coordinate(tr.eye - disp, &group_coord_low);
        split_coordinate(tr.eye + disp, &group_coord_high);

        for (int32_t zH = group_coord_low.z; zH <= group_coord_high.z; ++zH) {
        for (int32_t yH = group_coord_low.y; yH <= group_coord_high.y; ++yH) {
        for (int32_t xH = group_coord_low.x; xH <= group_coord_high.x; ++xH) {
            auto group_coord = glm::ivec3(xH, yH, zH);

            set_current_chunk_group(group_coord);
            if (current_chunk_group == nullptr) continue;

            float min_squared_dist;
            bool cull = cull_group(group_coord, &min_squared_dist);
            // Skip culled chunk groups or those so far away that they
            // cannot possibly contain chunks near enough to be drawn
            // using the mesh renderer.
            if (cull or min_squared_dist >= squared_thresh) continue;

            // Lookup the needed entry; scheduling if missing or dirty.
            // Draw the group if an entry for it is found.
            MeshEntry* entry = store.request_entry(group_coord);

            if (entry == nullptr) {
                store.enqueue(group_coord);
            }
            else {
                auto bit = renderer_mesh_dirty_flag;
                if (current_chunk_group->dirty_flags.load() & bit) {
                    store.enqueue(group_coord);
                    current_chunk_group->dirty_flags &= ~bit;
                }
                entries.emplace_back(entry, group_coord);
            }
        }
        }
        }

        draw_mesh_entries(entries);

        // Handle requests for new chunk groups.
        store.stage_from_queue(8);
        store.swap_in_from_staging(8);
    }

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
        const std::vector<std::pair<MeshEntry*, glm::ivec3>>& entries)
    {
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

    // Now time to write the AABB-raycasting renderer. Here goes...

    // Big-picture of data: I'm storing voxel data for chunk groups as
    // 3D textures. These textures live in RaycastEntry inside the
    // main cache of AsyncCache, below.
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
    // swap_in_from_staging) to give the compute shader time to
    // finish. See RaycastStaging::sync at time of writing.
    class RaycastStore : public AsyncCache<RaycastStaging, RaycastEntry>
    {
        // Back-pointer to the Renderer owning this RaycastStore.
        Renderer* owner = nullptr;

        static constexpr AsyncCacheArgs args = {
            4,          // modulus
            12,         // associativity (CHECK THIS IF THERE's FLICKERING)
            128,        // staging buffers
            2,          // worker threads
            10,         // condvar timeout in milliseconds
        };

        // OpenGL program for the earlier-mentioned compute shader.
        GLuint program_id = 0;

      public:
        RaycastStore(Renderer* owner_) :
            AsyncCache(args), owner(owner_)
        {
            program_id = make_program( { "staging.comp" } );
        }

        ~RaycastStore()
        {
            glDeleteProgram(program_id);
        }

        // Read from disk/memory data for the given chunk group and
        // fill the staging buffer's data with that of the chunk
        // group.
        bool stage(RaycastStaging* stage, glm::ivec3 group_coord) override
        {
            thread_local std::unique_ptr<ViewWorldCache> world_cache_ptr;

            if (world_cache_ptr == nullptr) {
                world_cache_ptr.reset(new ViewWorldCache(owner->world_handle));
            }

            auto& entry = world_cache_ptr->get_entry(group_coord);
            const BinChunkGroup* group_ptr = entry.chunk_group_ptr.get();
            if (group_ptr == nullptr) return false;

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
            return true;
        }

        // Swapping into the cache is done in two steps:
        //
        // 1. Dispatch the compute shader to transfer data from SSBO in
        // the staging buffer to the texture in the cache entry.
        //
        // 2. Wait for that compute shader to finish, then swap into
        // the cache.
        bool swap_in(
            RaycastStaging* staging,
            std::unique_ptr<RaycastEntry>* p_uptr_entry) override
        {
            // Step 1
            if (staging->sync == nullptr) {
                glUseProgram(program_id);
                // Hook up the SSBO binding to the SSBO index that will
                // be used to deliver voxel data.
                glShaderStorageBlockBinding(program_id,
                    chunk_group_voxels_program_index,
                    chunk_group_voxels_binding_index);
                // Hook up the out_image to the correct image unit.
                glUniform1i(staging_image_program_index, staging_image_unit);

                // Bind the correct texture image.
                if (*p_uptr_entry == nullptr) {
                    p_uptr_entry->reset(new RaycastEntry);
                }
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
    };

    // Render, to the current framebuffer, chunks around the camera
    // using the AABB-raycast algorithm. Chunks that are near
    // enough to have been drawn using the mesh algorithm will
    // not be re-drawn.
    void render_world_raycast_step() noexcept
    {
        RaycastStore& store = *raycast_store;
        store.begin_frame();

        // Resize the cache to a size suitable for the given render distance.
        auto far_plane = tr.far_plane;
        uint32_t modulus = uint32_t(std::max(4.0, ceil(far_plane / 128.0)));
        store.set_modulus(modulus, glFinish);

        // Count and limit the number of new chunk groups added to the
        // RaycastStore this frame.
        int remaining_new_chunk_groups_allowed =
            tr.max_frame_new_chunk_groups;

        // Collect all chunk groups needing to be raycast, and sort from
        // nearest to furthest (reverse painters). This fascilitates
        // the early depth test optimization.
        float squared_thresh = tr.far_plane;
        squared_thresh *= squared_thresh;
        std::vector<std::pair<float, glm::ivec3>> group_coord_by_depth;

        unsigned drawn_group_count = 0;
        glm::ivec3 group_coord_low, group_coord_high;
        glm::dvec3 disp(tr.far_plane);
        split_coordinate(tr.eye - disp, &group_coord_low);
        split_coordinate(tr.eye + disp, &group_coord_high);

        for (int32_t zH = group_coord_low.z; zH <= group_coord_high.z; ++zH) {
        for (int32_t yH = group_coord_low.y; yH <= group_coord_high.y; ++yH) {
        for (int32_t xH = group_coord_low.x; xH <= group_coord_high.x; ++xH) {
            auto group_coord = glm::ivec3(xH, yH, zH);
            set_current_chunk_group(group_coord);
            if (current_chunk_group == nullptr) continue;

            float min_squared_dist;
            bool cull = cull_group(group_coord, &min_squared_dist);
            if (cull or min_squared_dist >= squared_thresh) continue;
            group_coord_by_depth.emplace_back(min_squared_dist, group_coord);
            ++drawn_group_count;
        }
        }
        }

        if (!disable_zcull_sort) {
            auto lt_depth = [] (const auto& left, const auto& right)
            {
                return left.first < right.first;
            };
            std::sort(group_coord_by_depth.begin(),
                      group_coord_by_depth.end(),
                      lt_depth);
        }

        // Dispatch compute shaders mentioned earlier.
        raycast_store->swap_in_from_staging(80);

        // Now collect all the RaycastEntry objects we will draw,
        // scheduling loading new/dirty entries as we go along.
        std::vector<std::pair<RaycastEntry*, glm::ivec3>> entries;
        for (auto pair : group_coord_by_depth) {
            glm::ivec3 group_coord = pair.second;
            RaycastEntry* entry = raycast_store->request_entry(group_coord);

            if (entry == nullptr) {
                if (remaining_new_chunk_groups_allowed-- > 0) {
                    raycast_store->enqueue(group_coord);
                }
            }
            else {
                entries.emplace_back(entry, group_coord);

                set_current_chunk_group(group_coord);
                assert(current_chunk_group != nullptr);
                auto bit = renderer_raycast_dirty_flag;
                if ((current_chunk_group->dirty_flags.load() & bit)) {
                    store.enqueue(group_coord);
                    current_chunk_group->dirty_flags &= ~bit;
                }
            }
        }

        // Draw all the raycast entries.
        draw_raycast_entries(entries);

        // Now start filling the new staging buffers.
        PANIC_IF_GL_ERROR;
        raycast_store->stage_from_queue(80);
    }

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

    // Draw the given list of RaycastEntry/group coordinate pairs.
    void draw_raycast_entries(
        const std::vector<std::pair<RaycastEntry*, glm::ivec3>>& entries)
    {
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
        glUniform1i(u.chunk_debug, chunk_debug);
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

    void draw_frame()
    {
        tr = sync_camera_ptr->get_transforms_vk(); // See glClipControl
        int x = tr.frame_x, y = tr.frame_y;

        glViewport(0, 0, x, y);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
        if (tr.target_fragments > 0) {
            fprintf(stderr, "Not implemented: target_fragments.\n");
            // bind_global_f32_depth_framebuffer(x, y, tr.target_fragments);
        }
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
        render_world_mesh_step();
        render_world_raycast_step();
        render_background();
        if (tr.target_fragments > 0) {
            // finish_global_f32_depth_framebuffer(x, y);
        }

        window_ptr->swap_buffers();
    }

    static void render_loop(Renderer* renderer)
    {
        renderer->window_ptr->gl_make_current();

        renderer->mesh_store.reset(new MeshStore(renderer));
        renderer->raycast_store.reset(new RaycastStore(renderer));

        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glProvokingVertex(GL_FIRST_VERTEX_CONVENTION);
        glClipControl(GL_UPPER_LEFT, GL_ZERO_TO_ONE);
        glFrontFace(GL_CW);

        // In principle we just have to draw frames in a loop now, but
        // this function is bigger than you expect since I want to
        // compute FPS.

        // Seconds (since glfw initialization?) of previous drawn frame.
        double previous_update = glfwGetTime();

        double previous_fps_update = glfwGetTime();
        int frames = 0;
        double frame_time = 0;

        while (!renderer->thread_exit_flag.load()) {
            // First, actually draw the frame.
            renderer->draw_frame();

            // Calculate (approximate) frame time. Require at least 2
            // ms between frames (workaround to system freeze bug).
            double now;
            double dt;
            do {
                now = glfwGetTime();
                dt = now - previous_update;
            } while (dt < 0.002);
            previous_update = now;

            // Update FPS and frame time. Periodically report the
            // frame time back to the Renderer.
            ++frames;
            frame_time = std::max(frame_time, dt);
            if (now - previous_fps_update >= 0.5) {
                renderer->fps.store(frames / (now - previous_fps_update));
                renderer->frame_time.store(frame_time);

                previous_fps_update = now;
                frames = 0;
                frame_time = 0;
            }
        }
    }
};

Renderer* new_renderer(
    std::shared_ptr<Window> window,
    WorldHandle world,
    std::shared_ptr<SyncCamera> camera)
{
    return new Renderer(std::move(window), std::move(world), std::move(camera));
}

void delete_renderer(Renderer* renderer)
{
    delete renderer;
}

double get_fps(const Renderer& renderer)
{
    return renderer.fps.load();
}

double get_frame_time(const Renderer& renderer)
{
    return renderer.frame_time.load();
}

} // end namespace
