// Implementation of the hybrid voxel renderer. I don't see how I can
// make much of this exception safe, so I wrap everything in a
// noexcept at the end. As usual with OpenGL, none of this is
// thread-safe either.

#include "myricube.hh"

#include <algorithm>
#include <memory>
#include <stdio.h>
#include <typeinfo>
#include <utility>
#include <vector>

#include "camera.hh"
#include "chunk.hh"
#include "glad/glad.h"
#include "renderer.hh"
#include "shaders.hh"
#include "SDL2/SDL.h"

namespace myricube {

bool chunk_debug = false;

// My new plan for the GPU voxel mesh: each visible voxel will
// be represented by ONE vertex in a VBO. This vertex encodes,
// in a packed bitfield format,
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

    // Pack this vertex given the vertex's residue coordinates
    // (i.e. local coordinate in coordinate system of a chunk group)
    // and 24-bit color.
    //
    // This constructor defines the packed layout.
    MeshVoxelVertex(uint8_t x, uint8_t y, uint8_t z,
                    uint8_t red, uint8_t green, uint8_t blue,
                    uint32_t visible_face_bits)
    {
        static_assert(group_size <= 255,
                      "group too big for 8-bit unsigned coordinates.");
        assert(visible_face_bits == (visible_face_bits & all_face_bits));
        packed_residue_face_bits = uint32_t(x) << x_shift
                                 | uint32_t(y) << y_shift
                                 | uint32_t(z) << z_shift
                                 | visible_face_bits;
        packed_color = blue | uint32_t(green) << 8 | uint32_t(red) << 16;
    }
};

// The maximum number of MeshVoxelVerts needed for one chunk.  For now
// it's just the total number of voxels per chunk. I'm sure there's a
// lower possible bound but for now I'll be conservative (even though
// I'm desperate for GPU memory).
constexpr size_t chunk_max_verts = chunk_size * chunk_size * chunk_size;

// CPU-side copy of the mesh for a single chunk.
struct ChunkMesh
{
    MeshVoxelVertex verts[chunk_max_verts];
    size_t vert_count = 0;
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
    uint64_t world_id = 0;
    glm::ivec3 group_coord;

    // mesh_array[z][y][x] is the mesh for ChunkGroup::chunk_array[z][y][x].
    ChunkMesh mesh_array[edge_chunks][edge_chunks][edge_chunks];

    // OpenGL name for the VBO used to store the meshes of this chunk group.
    // 0 when not yet allocated.
    GLuint vbo_name = 0;

    // Size in bytes of the vbo's data store on the GPU.
    static constexpr GLsizeiptr vbo_bytes =
        GLsizeiptr(sizeof(ChunkMesh::verts)
                  * edge_chunks * edge_chunks * edge_chunks);

    // Return the offset (in number of MeshVoxelVertexes, not bytes)
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
        auto vert_sz = GLsizeiptr(sizeof(ChunkMesh::verts[0]));
        GLsizeiptr off = vert_offset(x, y, z) * vert_sz;
        assert(size_t(off) < vbo_bytes);
        return off;
    }

    MeshEntry() = default;

    // Get some RAII going.
    ~MeshEntry()
    {
        glDeleteBuffers(1, &vbo_name);
        vbo_name = 0;
        world_id = 0;
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

// All voxels within a chunk group share a single 3D texture on the
// GPU, as well as one VBO used to store the array of AABB for the
// chunks within the group. This is the handle for the GPU data needed
// to raycast one chunk group.
struct RaycastEntry
{
    // World ID and group coordinate of the chunk group this texture
    // and AABB array is for. This can be used to tell if this entry
    // is correct or needs to be replaced in RaycastStore.
    uint64_t world_id = 0;
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
        world_id = 0;
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
template <class Entry, uint32_t DefaultModulus, uint32_t Assoc>
class BaseStore
{
    struct CacheSet
    {
        uint64_t last_access[Assoc] = {0};
        Entry slots[Assoc];
    };
    // Monotonic increasing counter for least-recently accessed cache
    // eviction algorithm.
    uint64_t access_counter = 0;

    // Cache is (modulus x modulus x modulus x Assoc) in size.
    uint32_t modulus = DefaultModulus;

    // Conceptually, a 3D (modulus x modulus x modulus) array of cache
    // sets (CacheSet itself comprises the Assoc dimension).
    //
    // Declared after modulus defensively.
    std::unique_ptr<CacheSet[]> cache_sets {
        new CacheSet[DefaultModulus * DefaultModulus * DefaultModulus] };

    // Return a pointer to the location that the Entry for the chunk
    // group with the given group coordinate and world id should be.
    // The p_valid bool tells us whether the returned Entry matches
    // the one being searched for.
    Entry* cached_location(
        bool* p_valid, glm::ivec3 group_coord, uint64_t world_id)
    {
        uint32_t x = uint32_t(group_coord.x ^ 0x8000'0000) % modulus;
        uint32_t y = uint32_t(group_coord.y ^ 0x8000'0000) % modulus;
        uint32_t z = uint32_t(group_coord.z ^ 0x8000'0000) % modulus;
        CacheSet& cache_set = cache_sets[
            z * modulus * modulus + y * modulus + x];

        // Try to search for a valid Entry for the requested world &
        // group coordinate. Mark it as accessed if so
        for (unsigned i = 0; i < Assoc; ++i) {
            Entry* entry = &cache_set.slots[i];
            if (entry->world_id == world_id
            and entry->group_coord == group_coord) {
                cache_set.last_access[i] = ++access_counter;
                *p_valid = true;
                return entry;
            }
        }

        // Not found at this point; need to choose the least-recently
        // used entry to evict. This is much colder code than above.
        uint64_t min_access = cache_set.last_access[0];
        unsigned evict_idx = 0;
        for (unsigned i = 1; i < Assoc; ++i) {
            if (cache_set.last_access[i] < min_access) evict_idx = i;
        }
        Entry* entry = &cache_set.slots[evict_idx];

        *p_valid = false;
        cache_set.last_access[evict_idx] = ++access_counter;
        return entry;
    }

  public:
    uint64_t eviction_count = 0;
    // Return a valid, updated Entry for the given chunk group
    // extracted from the given world. TODO document templates needed
    // for this to work.
    //
    // If allow_failure is true, return null instead of evicting a
    // cache entry if the needed entry is not found.
    template <typename EntryFiller>
    Entry* update(PositionedChunkGroup& pcg,
                  VoxelWorld& world,
                  bool allow_failure=false)
    {
        glm::ivec3 gc = group_coord(pcg);
        bool valid;
        Entry* entry = cached_location(&valid, gc, world.id());

        // If the Entry in the array already corresponds to the given
        // chunk group; just update it (deal with dirty chunks) and
        // return.
        if (valid) {
            EntryFiller::update(pcg, world, entry);
            return entry;
        }

        if (allow_failure) return nullptr;

        EntryFiller::replace(pcg, world, entry);
        ++eviction_count;
        return entry;
    }

    // Shrink/grow the cache, moving every cache entry (belonging to
    // the specified world) to its correct location in the new cache.
    // Entries for other worlds are skipped over.
    //
    // I put a noexcept because I don't really know how to properly
    // handle exceptions in the middle of OpenGL code, not because
    // this is guaranteed foolproof.
    void set_modulus(uint32_t new_modulus, uint64_t world_id) noexcept
    {
        if (new_modulus == modulus) return;
        assert(int(new_modulus) > 0);

        const uint32_t old_modulus = modulus;
        modulus = new_modulus;

        // Swap the new and old caches.
        std::unique_ptr<CacheSet[]> old_cache_sets(
            new CacheSet[new_modulus * new_modulus * new_modulus]);
        swap(old_cache_sets, cache_sets);

        fprintf(stderr, "Set modulus to %i\n", int(modulus));

        // Put the old entries in their correct place.
        if (world_id != 0) {
            const size_t end = old_modulus * old_modulus * old_modulus;
            for (size_t i = 0; i < end; ++i) {
                CacheSet& old_set = old_cache_sets[i];
                for (size_t a = 0; a < Assoc; ++a) {
                    Entry& old_entry = old_set.slots[a];
                    if (old_entry.world_id != world_id) continue;

                    bool valid;
                    Entry* new_location = cached_location(
                        &valid, old_entry.group_coord, world_id);
                    swap(*new_location, old_entry);
                    // Note: May want to preserve the old last_accessed time.
                    // For now, each cached_location updates the access time,
                    // so earlier repositioned cache entries are arbitrarily
                    // considered lower priority later.
                }
            }
        }
    }
};

class MeshStore : public BaseStore<MeshEntry, 2, 8> { };
class RaycastStore : public BaseStore<RaycastEntry, 3, 12> { };

// AABB is drawn as a unit cube from (0,0,0) to (1,1,1), which is
// stretched and positioned to the right shape and position in space.
// The normals are "flat shaded" and needed for dealing with
// floor/ceil rounding errors. (Basically, I add a bit of the normal
// vector to each face in order to ensure the voxels fit comfortably
// within the AABB).
//
// [positions] [normal]
//
// This really should be an array of structs, but I copied this from
// some older project and I'm too scared to screw it up.
static const float unit_box_vertices[48] =
{
    0, 1, 1,   -1, 0, 0,    // 0 -x face provoking vertex
    0, 0, 1,   0, 0, 1,     // 1 +z face provoking vertex
    1, 0, 1,   1, 0, 0,     // 2 +x face provoking vertex
    1, 1, 1,   0, 0, 0,     // 3 unused as provoking vertex
    0, 1, 0,   0, 1, 0,     // 4 +y face provoking vertex
    0, 0, 0,   0, -1, 0,    // 5 -y face provoking vertex
    1, 0, 0,   0, 0, 0,     // 6 unused as provoking vertex
    1, 1, 0,   0, 0, -1,    // 7 -z face provoking vertex
};

static const GLushort unit_box_elements[36] = {
    6, 7, 2, 7, 3, 2,   // +x face
    4, 5, 0, 5, 1, 0,   // -x face
    0, 3, 4, 3, 7, 4,   // +y face
    6, 2, 5, 2, 1, 5,   // -y face
    2, 3, 1, 3, 0, 1,   // +z face
    6, 5, 7, 5, 4, 7,   // -z face
};


// Voxels rendered using the mesh renderer also need a unit cube, but
// this time, instead of a normal I need a bit indicating which face
// each vertex is part of (so I can discard hidden faces).
struct VoxelUnitBoxVertex
{
    int32_t face_bit;
    float x, y, z;
};

static const VoxelUnitBoxVertex voxel_unit_box_vertices[24] =
{
    { neg_x_face_bit, 0, 1, 1 },
    { neg_x_face_bit, 0, 1, 0 },
    { neg_x_face_bit, 0, 0, 1 },
    { neg_x_face_bit, 0, 0, 0 },

    { pos_x_face_bit, 1, 0, 0 },
    { pos_x_face_bit, 1, 1, 0 },
    { pos_x_face_bit, 1, 0, 1 },
    { pos_x_face_bit, 1, 1, 1 },

    { neg_y_face_bit, 0, 0, 0 },
    { neg_y_face_bit, 1, 0, 1 },
    { neg_y_face_bit, 0, 0, 1 },
    { neg_y_face_bit, 1, 0, 0 },

    { pos_y_face_bit, 1, 1, 1 },
    { pos_y_face_bit, 1, 1, 0 },
    { pos_y_face_bit, 0, 1, 0 },
    { pos_y_face_bit, 0, 1, 1 },

    { neg_z_face_bit, 0, 0, 0 },
    { neg_z_face_bit, 0, 1, 0 },
    { neg_z_face_bit, 1, 0, 0 },
    { neg_z_face_bit, 1, 1, 0 },

    { pos_z_face_bit, 0, 1, 1 },
    { pos_z_face_bit, 1, 0, 1 },
    { pos_z_face_bit, 1, 1, 1 },
    { pos_z_face_bit, 0, 0, 1 },
};

static const GLushort voxel_unit_box_elements[36] =
{
    0, 1, 2,    1, 3, 2,
    4, 5, 6,    5, 7, 6,
    8, 9, 10,   8, 11, 9,
    12, 13, 14, 15, 12, 14,
    16, 17, 18, 19, 18, 17,
    20, 21, 22, 21, 20, 23,
};

static bool culling_freeze = false;
static glm::mat4 frozen_vp;
static glm::ivec3 frozen_eye_group;
static glm::vec3 frozen_eye_residue;

// Renderer class. Instantiate it with the camera and world to render
// and use it once.
class Renderer
{
    Camera& camera;
    VoxelWorld& world;
  public:
    Renderer(VoxelWorld& world_, Camera& camera_) :
        camera(camera_), world(world_) { }

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
    bool cull_group(PositionedChunkGroup& pcg, float* squared_dist=nullptr)
    {
        if (group(pcg).total_visible == 0) return true;

        glm::mat4 vp = culling_freeze ? frozen_vp : camera.get_residue_vp();
        glm::ivec3 eye_group;
        glm::vec3 eye_residue;
        if (culling_freeze) {
            eye_group = frozen_eye_group;
            eye_residue = frozen_eye_residue;
        } else {
            camera.get_eye(&eye_group, &eye_residue);
        }

        // Never cull the group the eye is in (this avoids
        // pathological cases e.g. 7 corners behind eye and 1 in front
        // and out of the frustum).
        if (group_coord(pcg) == eye_group) {
            if (squared_dist) *squared_dist = 0.0f;
            return false;
        }

        eye_residue = glm::floor(eye_residue); // to match decide_chunk.

        // Position of this chunk group relative to the group that the
        // eye is in, in voxel units.
        glm::vec3 low_corner(group_size * (group_coord(pcg) - eye_group));

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
                     Chunk& chunk)
    {
        if (chunk.total_visible == 0) return cull;
        glm::ivec3 aabb_low, aabb_high;
        chunk.get_aabb(&aabb_low, &aabb_high);
        glm::ivec3 eye_group;
        glm::vec3 eye_residue;
        camera.get_eye(&eye_group, &eye_residue);
        auto far_plane = camera.get_far_plane();
        auto raycast_thresh = camera.get_raycast_threshold();

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

    // Given a chunk, fill in the mesh needed to render it. I also
    // need the residue coordinate of the lower-left of the chunk
    // (i.e.  the position of the chunk relative to the chunk group it
    // is in) because the mesh is defined to be in residue
    // coordinates.
    static void fill_mesh_verts(ChunkMesh& mesh,
                                const Chunk& chunk,
                                glm::ivec3 chunk_residue)
    {
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
            return chunk(coord).visible;
        };

        auto visit_voxel = [&mesh, &chunk, visible_block, chunk_residue]
        (glm::ivec3 coord)
        {
            Voxel v = chunk(coord);
            if (!v.visible) return;

            auto r = v.red;
            auto g = v.green;
            auto b = v.blue;
            uint8_t x = uint8_t(coord.x) + chunk_residue.x;
            uint8_t y = uint8_t(coord.y) + chunk_residue.y;
            uint8_t z = uint8_t(coord.z) + chunk_residue.z;

            MeshVoxelVertex vert(x, y, z, r, g, b, 0);

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
                assert(mesh.vert_count < chunk_max_verts);
                mesh.verts[mesh.vert_count++] = vert;
            }
        };

        mesh.vert_count = 0;

        for (int z = 0; z < chunk_size; ++z) {
            for (int y = 0; y < chunk_size; ++y) {
                for (int x = 0; x < chunk_size; ++x) {
                    visit_voxel(glm::ivec3(x, y, z));
                }
            }
        }
    }

    // Fill the given MeshEntry with mesh data for the given chunk
    // group from the given world.
    static void replace(PositionedChunkGroup& pcg,
                        VoxelWorld& world,
                        MeshEntry* entry)
    {
        entry->world_id = world.id();
        entry->group_coord = group_coord(pcg);

        if (entry->vbo_name == 0) {
            glGenBuffers(1, &entry->vbo_name);
            glBindBuffer(GL_ARRAY_BUFFER, entry->vbo_name);
            glBufferStorage(GL_ARRAY_BUFFER,
                            MeshEntry::vbo_bytes,
                            nullptr,
                            GL_DYNAMIC_STORAGE_BIT);
            PANIC_IF_GL_ERROR;
        }

        // Parasitically depend on the update, using the AlwaysDirty
        // parameter to get every chunk's mesh copied onto the new (or
        // reused) VBO even if not marked dirty.
        update<true>(pcg, world, entry);
    }

    // Update the given mesh entry with new data from dirty chunks.
    //
    // Requires that the MeshEntry corresponds to the given chunk
    // group of the given world.
    template <bool AlwaysDirty=false>
    static void update(PositionedChunkGroup& pcg,
                       VoxelWorld& world,
                       MeshEntry* entry)
    {
        assert(entry->world_id == world.id());
        assert(entry->group_coord == group_coord(pcg));

        auto vbo_name = entry->vbo_name;
        assert(vbo_name != 0);

        // Recompute and reupload the mesh for any dirty chunk, and
        // reset the chunk dirty flag.
        for (int z = 0; z < edge_chunks; ++z) {
            for (int y = 0; y < edge_chunks; ++y) {
                for (int x = 0; x < edge_chunks; ++x) {
                    Chunk& chunk = group(pcg).chunk_array[z][y][x];
                    if (!AlwaysDirty and !chunk.mesh_dirty) continue;
                    if (chunk.total_visible == 0) continue;

                    ChunkMesh& mesh = entry->mesh_array[z][y][x];
                    glm::ivec3 residue = glm::ivec3(x,y,z) * chunk_size;
                    fill_mesh_verts(mesh, chunk, residue);
                    glNamedBufferSubData(vbo_name,
                                         entry->byte_offset(x, y, z),
                                         sizeof mesh.verts[0] * mesh.vert_count,
                                         mesh.verts);
                    chunk.mesh_dirty = false;
                }
            }
        }
        PANIC_IF_GL_ERROR;
    }

    // Render, to the current framebuffer, chunks near the camera
    // using the conventional mesh-based algorithm.
    //
    // Each MeshEntry contains a list of visible voxels: their residue
    // coordinates, colors, and which of their 6 faces are visible.
    // My plan is to use instanced rendering to draw a bunch of unit
    // cubes in the correct locations, using degenerate triangles to
    // hide the hidden faces.
    void render_world_mesh_step() noexcept
    {
        MeshStore& store = camera.get_mesh_store();
        store.eviction_count = 0;

        glm::mat4 residue_vp_matrix = camera.get_residue_vp();
        glm::vec3 eye_residue;
        glm::ivec3 eye_group;
        camera.get_eye(&eye_group, &eye_residue);

        static GLuint vao = 0;
        static GLuint program_id;
        static GLint mvp_matrix_idx;

        // Position of camera eye relative to the origin of the group
        // (origin == group_size times the group coordinate).
        static GLint eye_relative_group_origin_id;

        static GLint far_plane_squared_id;
        static GLint fog_enabled_id;

        if (vao == 0) {
            program_id = make_program({ "mesh.vert", "mesh.frag" });
            mvp_matrix_idx = glGetUniformLocation(program_id, "mvp_matrix");
            assert(mvp_matrix_idx >= 0);
            eye_relative_group_origin_id = glGetUniformLocation(program_id,
                "eye_relative_group_origin");
            assert(eye_relative_group_origin_id >= 0);
            far_plane_squared_id = glGetUniformLocation(program_id,
                "far_plane_squared");
            assert(far_plane_squared_id >= 0);
            fog_enabled_id = glGetUniformLocation(program_id,
                "fog_enabled");
            assert(fog_enabled_id >= 0);
            glGenVertexArrays(1, &vao);
            glBindVertexArray(vao);
            PANIC_IF_GL_ERROR;

            // The above was all boilerplate crap; this is where I
            // bind once and forever the unit box buffers to the VAO.
            GLuint buffers[2];
            glGenBuffers(2, buffers);
            glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[1]);

            glBufferStorage(GL_ARRAY_BUFFER,
                            sizeof voxel_unit_box_vertices,
                            voxel_unit_box_vertices,
                            0);
            glBufferStorage(GL_ELEMENT_ARRAY_BUFFER,
                            sizeof voxel_unit_box_elements,
                            voxel_unit_box_elements,
                            0);

            glVertexAttribPointer(
                unit_box_vertex_idx,
                3,
                GL_FLOAT,
                false,
                sizeof(VoxelUnitBoxVertex),
                (void*) offsetof(VoxelUnitBoxVertex, x));
            glEnableVertexAttribArray(unit_box_vertex_idx);
            glVertexAttribIPointer(
                unit_box_face_bit_idx,
                1,
                GL_UNSIGNED_INT,
                sizeof(VoxelUnitBoxVertex),
                (void*) offsetof(VoxelUnitBoxVertex, face_bit));
            glEnableVertexAttribArray(unit_box_face_bit_idx);
            PANIC_IF_GL_ERROR;
        }

        // If we're not using D-tier Intel drivers, this should
        // restore the element array binding above...
        glBindVertexArray(vao);

        glUseProgram(program_id);
        auto far_plane = camera.get_far_plane();
        glUniform1i(far_plane_squared_id, far_plane * far_plane);
        glUniform1i(fog_enabled_id, camera.get_fog());
        PANIC_IF_GL_ERROR;

        auto draw_group = [&] (PositionedChunkGroup& pcg)
        {
            if (group(pcg).total_visible == 0) return;
            MeshEntry* entry = store.update<Renderer>(pcg, world);

            assert(entry->vbo_name != 0);
            glBindBuffer(GL_ARRAY_BUFFER, entry->vbo_name);

            // Bind instanced (per-voxel) attributes.
            // TODO: Maybe make one VAO per MeshEntry?
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
            glm::vec3 model_offset = glm::vec3(group_coord(pcg) - eye_group)
                                   * float(group_size);
            glm::mat4 m = glm::translate(glm::mat4(1.0f), model_offset);
            glm::mat4 mvp = residue_vp_matrix * m;
            glUniformMatrix4fv(mvp_matrix_idx, 1, 0, &mvp[0][0]);

            // Similarly, the eye residue needs to be shifted by the
            // group's position.
            glm::vec3 eye_relative_group_origin = eye_residue - model_offset;
            glUniform3fv(eye_relative_group_origin_id, 1,
                &eye_relative_group_origin[0]);

            auto gc = group_coord(pcg);
            for (int z = 0; z < edge_chunks; ++z) {
                for (int y = 0; y < edge_chunks; ++y) {
                    for (int x = 0; x < edge_chunks; ++x) {
                        Chunk& chunk = group(pcg).chunk_array[z][y][x];
                        if (decide_chunk(gc, chunk) != draw_mesh) continue;

                        const ChunkMesh& mesh = entry->mesh_array[z][y][x];
                        if (mesh.vert_count == 0) continue;

                        // Longest function name EVER. I need it because
                        // every chunk's data is at a different offset
                        // within the VBO for instanced data (position/color).
                        glDrawElementsInstancedBaseVertexBaseInstance(
                            GL_TRIANGLES,
                            36, // 12 triangles, 3 verts each.
                            GL_UNSIGNED_SHORT,
                            nullptr,
                            mesh.vert_count,
                            0, // base vertex is always 0, for the unit box.
                            entry->vert_offset(x, y, z)
                            // ^^^ Instance offset, depends on chunk.
                        );
                    }
                }
            }
            PANIC_IF_GL_ERROR;
        };

        float squared_thresh = camera.get_raycast_threshold();
        squared_thresh *= squared_thresh;

        // TODO: Avoid visiting too-far-away chunk groups.
        for (PositionedChunkGroup& pcg : world.group_map) {
            float min_squared_dist;
            bool cull = cull_group(pcg, &min_squared_dist);
            // Skip culled chunk groups or those so far away that they
            // cannot possibly contain chunks near enough to be drawn
            // using the mesh renderer.
            if (cull or min_squared_dist >= squared_thresh) continue;
            draw_group(pcg);
        }

        if (store.eviction_count >= 5) {
            fprintf(stderr, "%i MeshStore evictions.\n",
                int(store.eviction_count));
        }

        glBindVertexArray(0);
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
                    bool aabb_dirty = always_dirty | chunk.aabb_gpu_dirty;
                    vbo_dirty |= aabb_dirty;

                    if (aabb_dirty) {
                        glm::ivec3 aabb_low, aabb_high;
                        chunk.get_aabb(&aabb_low, &aabb_high);
                        chunk.aabb_gpu_dirty = false;
                        entry->aabb_array[z][y][x] =
                            PackedAABB(aabb_low, aabb_high);
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
    void render_world_raycast_step() noexcept
    {
        RaycastStore& store = camera.get_raycast_store();
        store.eviction_count = 0;

        // Resize the cache to a size suitable for the given render distance.
        auto far_plane = camera.get_far_plane();
        uint32_t modulus = uint32_t(std::max(3.0, ceil(far_plane / 128.0)));
        store.set_modulus(modulus, world.id());

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
        // (origin == group_size times the group coordinate).
        static GLint eye_relative_group_origin_id;
        static GLint far_plane_squared_id;
        static GLint raycast_thresh_squared_id;
        static GLint fog_enabled_id;
        static GLint chunk_blocks_id;
        static GLint chunk_debug_id;

        if (vao == 0) {
            program_id = make_program({ "raycast.vert", "raycast.frag" });
            mvp_matrix_id = glGetUniformLocation(program_id, "mvp_matrix");
            assert(mvp_matrix_id >= 0);
            eye_relative_group_origin_id = glGetUniformLocation(program_id,
                "eye_relative_group_origin");
            assert(eye_relative_group_origin_id >= 0);
            chunk_blocks_id = glGetUniformLocation(program_id, "chunk_blocks");
            assert(chunk_blocks_id >= 0);
            chunk_debug_id = glGetUniformLocation(program_id, "chunk_debug");
            assert(chunk_debug_id >= 0);
            far_plane_squared_id = glGetUniformLocation(program_id,
                "far_plane_squared");
            assert(far_plane_squared_id >= 0);
            raycast_thresh_squared_id = glGetUniformLocation(program_id,
                "raycast_thresh_squared");
            assert(raycast_thresh_squared_id >= 0);
            fog_enabled_id = glGetUniformLocation(program_id,
                "fog_enabled");
            assert(fog_enabled_id >= 0);

            glGenVertexArrays(1, &vao);
            glBindVertexArray(vao);

            // The binding points for the non-instanced unit cube
            // vertices, and the element buffer binding, will be
            // stored forever in the VAO.
            glGenBuffers(1, &vertex_buffer_id);
            glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_id);
            glBufferStorage(
                GL_ARRAY_BUFFER, sizeof unit_box_vertices,
                unit_box_vertices, 0);

            glGenBuffers(1, &element_buffer_id);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer_id);
            glBufferStorage(
                GL_ELEMENT_ARRAY_BUFFER, sizeof unit_box_elements,
                unit_box_elements, 0);

            glVertexAttribPointer(
                unit_box_vertex_idx,
                3,
                GL_FLOAT,
                false,
                sizeof(float) * 6,
                (void*)0);
            glEnableVertexAttribArray(unit_box_vertex_idx);

            glVertexAttribPointer(
                unit_box_normal_idx,
                3,
                GL_FLOAT,
                false,
                sizeof(float) * 6,
                (void*)12);
            glEnableVertexAttribArray(unit_box_normal_idx);
        }

        glBindVertexArray(vao);
        glUseProgram(program_id);
        glUniform1i(chunk_debug_id, chunk_debug);
        glUniform1i(fog_enabled_id, camera.get_fog());
        glActiveTexture(GL_TEXTURE0);
        glUniform1i(chunk_blocks_id, 0);
        glUniform1i(far_plane_squared_id, far_plane * far_plane);
        auto raycast_thr = camera.get_raycast_threshold();
        glUniform1i(raycast_thresh_squared_id, raycast_thr * raycast_thr);
        PANIC_IF_GL_ERROR;

        std::vector<PositionedChunkGroup*> deferred_pcg;

        // My plan is to use instanced rendering to draw the chunk AABBs
        // of this chunk group. The base box is a 1x1x1 unit box, which
        // is stretched and repositioned in the vertex shader to the
        // true AABB.
        //
        // This is a bit of an experimental after-the-fact hack but
        // the allow_defer and deferred_pcg variables are used to
        //
        // A. Throttle the number of RaycastStore updates done subject
        //    to the limits of camera.max_raycast_evict.
        //
        // B. Put some time between updating a missing RaycastEntry and
        //    drawing its corresponding chunk group. This /in principle/
        //    prevents texture updates from stalling the GPU so much.
        auto draw_group = [&] (PositionedChunkGroup& pcg, bool allow_defer)
        {
            // if (group(pcg).total_visible == 0) return;
            RaycastEntry* entry = store.update<Renderer>(
                pcg, world, allow_defer);

            if (entry == nullptr) {
                if (deferred_pcg.size() < camera.max_raycast_evict) {
                    store.update<Renderer>(pcg, world, false);
                    deferred_pcg.push_back(&pcg);
                }
                return;
            }
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

            // The unit box vertex attribs should already be bound by VAO.
            //
            // Get the instanced vertex attribs going (i.e. bind the
            // AABB residue coords, packed as integers). I should
            // probably make a VAO per chunk group in the
            // RaycastEntry.
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

        // Collect all chunk groups needing to be raycast, and sort from
        // nearest to furthest (reverse painters). This fascilitates
        // the early depth test optimization.
        float squared_thresh = camera.get_far_plane();
        squared_thresh *= squared_thresh;
        std::vector<std::pair<float, PositionedChunkGroup*>> pcg_by_depth;

        // TODO: Avoid visiting far away chunk groups.
        for (PositionedChunkGroup& pcg : world.group_map) {
            float min_squared_dist;
            bool cull = cull_group(pcg, &min_squared_dist);
            if (cull or min_squared_dist >= squared_thresh) continue;
            pcg_by_depth.emplace_back(min_squared_dist, &pcg);
        }

        std::sort(pcg_by_depth.begin(), pcg_by_depth.end());
        for (auto& pair : pcg_by_depth) {
            draw_group(*pair.second, true);
        }
        for (PositionedChunkGroup* p_pcg : deferred_pcg) {
            draw_group(*p_pcg, false);
        }

        if (store.eviction_count >= 5) {
            fprintf(stderr, "%i RaycastStore evictions.\n",
                int(store.eviction_count));
        }

        glBindVertexArray(0);
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

void toggle_culling_freeze(Camera& current_camera)
{
    if (culling_freeze) {
        culling_freeze = false;
    }
    else {
        culling_freeze = true;
        frozen_vp = current_camera.get_residue_vp();
        current_camera.get_eye(&frozen_eye_group, &frozen_eye_residue);
    }
}

static GLuint f32_depth_framebuffer = 0;
static GLuint f32_depth_texture = 0;
static GLuint f32_depth_renderbuffer = 0;
static int f32_depth_framebuffer_x, f32_depth_framebuffer_y = 0;

void bind_global_f32_depth_framebuffer(int screen_x, int screen_y)
{
    // Destroy the framebuffer and its attachments if the screen size changed.
    if (screen_x != f32_depth_framebuffer_x
    or screen_y != f32_depth_framebuffer_y) {
        if (f32_depth_framebuffer != 0) {
            fprintf(stderr, "Destroying old %i x %i framebuffer.\n",
                f32_depth_framebuffer_x, f32_depth_framebuffer_y);
            glDeleteFramebuffers(1, &f32_depth_framebuffer);
            glDeleteTextures(1, &f32_depth_texture);
            glDeleteRenderbuffers(1, &f32_depth_renderbuffer);
        }

        f32_depth_framebuffer = 0;
        f32_depth_framebuffer_x = screen_x;
        f32_depth_framebuffer_y = screen_y;
        PANIC_IF_GL_ERROR;
    }

    // (Re-) create the framebuffer if needed, and bind it.
    if (f32_depth_framebuffer == 0) {
        // Create framebuffer.
        fprintf(stderr, "Creating %i x %i framebuffer.\n",
                f32_depth_framebuffer_x, f32_depth_framebuffer_y);
        glCreateFramebuffers(1, &f32_depth_framebuffer);
        glBindFramebuffer(GL_FRAMEBUFFER, f32_depth_framebuffer);
        PANIC_IF_GL_ERROR;

        // Add depth buffer.
        glCreateRenderbuffers(1, &f32_depth_renderbuffer);
        glBindRenderbuffer(GL_RENDERBUFFER, f32_depth_renderbuffer); // XXX
        glNamedRenderbufferStorage(f32_depth_renderbuffer,
                                   GL_DEPTH_COMPONENT, // XXX
                                   f32_depth_framebuffer_x,
                                   f32_depth_framebuffer_y);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER,
                                  GL_DEPTH_ATTACHMENT,
                                  GL_RENDERBUFFER,
                                  f32_depth_renderbuffer);
        PANIC_IF_GL_ERROR;

        // Add color buffer.
        glCreateTextures(GL_TEXTURE_2D, 1, &f32_depth_texture);
        glBindTexture(GL_TEXTURE_2D, f32_depth_texture); // XXX
        glTextureParameteri(f32_depth_texture, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTextureParameteri(f32_depth_texture, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTextureParameteri(f32_depth_texture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTextureParameteri(f32_depth_texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTextureStorage2D(f32_depth_texture, 1, GL_RGBA8,
                           f32_depth_framebuffer_x, f32_depth_framebuffer_y);
        glFramebufferTexture2D(GL_FRAMEBUFFER,
                               GL_COLOR_ATTACHMENT0,
                               GL_TEXTURE_2D,
                               f32_depth_texture,
                               0);
        PANIC_IF_GL_ERROR;

        // Wire up the only color output.
        GLenum tmp = GL_COLOR_ATTACHMENT0;
        glDrawBuffers(1, &tmp);

        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        {
            panic("Incomplete 32-bit depth framebuffer: "
                 + std::to_string(glCheckFramebufferStatus(GL_FRAMEBUFFER)));
        }
    }
    else {
        glBindFramebuffer(GL_FRAMEBUFFER, f32_depth_framebuffer);
    }
    PANIC_IF_GL_ERROR;
}

void finish_global_f32_depth_framebuffer()
{
    assert(f32_depth_framebuffer != 0);
    glBlitNamedFramebuffer(
        f32_depth_framebuffer,
        0, // Write to window framebuffer
        0, 0,
        f32_depth_framebuffer_x, f32_depth_framebuffer_y,
        0, 0,
        f32_depth_framebuffer_x, f32_depth_framebuffer_y,
        GL_COLOR_BUFFER_BIT,
        GL_NEAREST);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glViewport(0, 0, f32_depth_framebuffer_x, f32_depth_framebuffer_y);
    PANIC_IF_GL_ERROR;
}

void enable_debug_callback()
{
    fprintf(stderr, "Registering debug callback.\n");
    auto callback = [] (GLenum source, GLenum type, GLuint id,
                        GLenum severity, GLsizei length,
                        const GLchar* message, const void* userParam)
    {
        fputs(message, stderr);
    };

    glDebugMessageCallback(callback, nullptr);
}

} // end namespace
