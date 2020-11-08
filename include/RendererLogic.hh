// Abstract base class for implementing a renderer thread.  Basically,
// this is all the "business logic" for rendering (i.e. everything
// that doesn't have to touch the actual graphics API used). Inherit
// from this and implement the abstract functions using the actual
// graphics API to complete the renderer.
//
// Actual use of this class should be through the RenderThread class,
// defined elsewhere.
//
// High-level description: The renderer is attached to a window,
// camera, and world upon instantiation [see RenderArgs], and runs on
// its own thread in a loop [RendererLogic::render_loop]. The
// instantiation is performed in the same thread that actually runs
// the loop, so you can actually use constructors to do useful work.
//
// The main purpose of this class is to manage the GPU caches of chunk
// groups (see AsyncCache), load chunk groups from disk/memory and add
// them to the caches as-needed, and extract lists of chunk groups
// from the caches (culled based on camera location: cull_group and
// decide_chunk) and pass them on to the subclass's functions to
// render them.
//
// By "manage", I mean "delegate all actual API work to the subclass's
// functions" (look for the protected section to see their
// specifics). Also, when I say "lists of chunk group", what I really
// mean is lists of the MeshEntry/RaycastEntry template parameters,
// which should hold the 3D-textures or whatever needed to actually
// render a chunk group.

#ifndef MYRICUBE_RENDERERLOGIC_HH_
#define MYRICUBE_RENDERERLOGIC_HH_

#include "myricube.hh"

#include <algorithm>
#include <cassert>
#include <memory>
#include <thread>
#include <vector>

#include "AsyncCache.hh"
#include "camera.hh"
#include "MeshVoxelVertex.hh"
#include "PackedAABB.hh"
#include "RenderThread.hh"
#include "voxels.hh"

namespace myricube {

// (roughly) minimum distance from the camera that a chunk needs
// to be to switch from mesh to raycast graphics.
// Keep as int to avoid rounding errors in distance culling.
// Might make this configurable some day.
constexpr int raycast_threshold = 160;

// Slightly higher threshold than raycast_threshold. Chunk groups
// within this distance to the camera have their meshes loaded even
// if no part is actually drawn as a mesh (this is needed now that
// MeshStore has some latency in loading meshes).
constexpr int mesh_load_threshold = 20 + raycast_threshold;

// See the protected virtual section to see the functions you have to
// implement.
//
// MeshEntry: a single chunk group's data that is passed to the mesh
// rendering function (e.g. this could be a list of triangles)
//
// RaycastEntry: a single chunk group's data that is passed to the
// raycast rendering function (e.g. this could be a 3D texture).
//
// Both are stored in a 3D cache, which could be read by rendering
// functions. To avoid modifying them while they're used for
// rendering, there's also "staging" versions of the two above types.
// Data flows from a CPU-side chunk group to a staging buffer on a
// worker thread (see worker_stage), then from a staging buffer to the
// main cache on the main thread (see swap_in). It is up to you how to
// divide work between the workers and main thread.
//
// More specifically, MeshEntry and MeshShading are the EntryT and
// StagingT type parameters for the AsyncCache used for mesh rendering
// (i.e. MeshEntry should be the API data needed to actually render a
// chunk group with the implemented mesh renderer). RaycastEntry and
// RaycastStaging are EntryT and StagingT for raycast rendering.
//
// Throughout, "main thread" refers to the thread that created this
// RendererLogic and is running its render_loop.
template <typename MeshEntry,
          typename MeshStaging,
          typename RaycastEntry,
          typename RaycastStaging>
class RendererLogic : public RendererBase
{
  protected:
    // Window, camera, world arguments used to launch the render thread.
    const std::shared_ptr<Window> p_window;
    const std::shared_ptr<SyncCamera> p_camera;
    const WorldHandle world_handle;

    // Updated from SyncCamera at the start of every frame.
    CameraTransforms transforms;

  private:
    // Implementations of AsyncCache, manages device storage for
    // storing chunk groups. These are not initialized until after the
    // subclass constructor finishes (so that needed resources have a
    // chance to be initialized).
    //
    // NOTE: A special hack is done at RendererBase::destroy_stores
    // ensure these are destroyed BEFORE the derived class destructor
    // is run (but after the last frame).  This is needed in case
    // MeshEntry/RaycastEntry depends on the context set up by the
    // derived class.
    class MeshStore;
    class RaycastStore;
    std::unique_ptr<MeshStore> mesh_store;
    std::unique_ptr<RaycastStore> raycast_store;

    // Constructed from world_handle to efficiently view chunk groups
    // of the world. Access only through set_current_chunk_group.
    ViewWorldCache world_cache;

    // A bit hacky, but this is the last chunk group requested from
    // world_cache. Needed because new accesses invalidate
    // previously-returned pointers from world_cache.
    const BinChunkGroup* current_chunk_group = nullptr;

  protected:
    // Base class constructor really just copies stuff in.
    RendererLogic(RenderThread* p_back_, RenderArgs args) :
        RendererBase(p_back_),
        p_window(args.p_window),
        p_camera(args.p_camera),
        world_handle(args.world_handle),
        mesh_store(nullptr),
        raycast_store(nullptr),
        world_cache(args.world_handle),
        current_chunk_group(nullptr)
    {

    }

    virtual ~RendererLogic() = default;

    // These are the functions you need to implement.

    // Called by the main thread at the start of a new frame.
    virtual void begin_frame() = 0;

    // Called by the main thread. Given a vector of chunk groups
    // (represented as pairs of their MeshEntry and group coordinate),
    // draw them using the mesh rendering method.
    virtual void
    draw_mesh_entries(
        const std::vector<std::pair<MeshEntry*, glm::ivec3>>&) = 0;

    // Called by the main thread after draw_mesh_entries.
    //
    // Given a vector of chunk groups (represented as pairs of their
    // RaycastEntry and group coordinate), draw them using the raycast
    // rendering method.
    //
    // NOTE: For things to look right, you'll have to use decide_chunk
    // to only draw chunks within the given chunk groups that are
    // within raycasting range.
    virtual void
    draw_raycast_entries(
        const std::vector<std::pair<RaycastEntry*, glm::ivec3>>&) = 0;

    // Called by the main thread at the end of a frame.
    virtual void end_frame() = 0;

    // Called by a worker thread. Read the given chunk group and
    // schedule whatever work is needed to start filling the given
    // MeshStaging with data needed to render the given chunk group
    // with the mesh rendering algorithm.
    //
    // WARNING: the BinChunkGroup pointer will become invalid as soon
    // as this function returns (be careful with asynchronous
    // transfers), and you need not clear the dirty flag of the chunk
    // group. XXX TODO make sure I implement this correctly.
    virtual void worker_stage(MeshStaging*, const BinChunkGroup*) = 0;

    // Like the above, but for raycasting. Same warning applies.
    virtual void worker_stage(RaycastStaging*, const BinChunkGroup*) = 0;

    // NOTE: Although the staging buffers (MeshStaging[] and
    // RaycastStaging[]) are filled by worker threads, they are
    // constructed and destroyed by the main thread.

    // Called by the main thread. Attempt to fill the given MeshEntry
    // with data staged in MeshStaging. (The MeshEntry is passed as a
    // pointer-to-unique-ptr so that you can exchange the MeshEntry
    // for an entirely new MeshEntry if needed).
    //
    // If you return true, MeshStaging MUST be ready for re-use as
    // soon as this function returns.
    //
    // Return true if successful, false if we need to retry later.
    virtual bool swap_in(MeshStaging*, std::unique_ptr<MeshEntry>*) = 0;

    // Same, but for raycasting.
    virtual bool swap_in(RaycastStaging*, std::unique_ptr<RaycastEntry>*) = 0;

    // NOTE: I really hope the above stays up-to-date if I change the
    // AsyncCache...

    // A bit of a wart, just pass in glFinish or vkQueueWaitIdle.
    // This is needed when we resize the cache due to render distance
    // change (note: only the main 3D cache, which is accessed by the
    // main thread only, is affected; the staging buffers shared with
    // the workers are not affected).
    virtual void main_thread_wait_idle() { wait_idle(); }

    // Like the above, but called after the worker threads are stopped
    // but before destroying the cache. (Difference from above is
    // destruction affects the staging buffers as well).
    virtual void wait_idle() = 0;

  private:
    // Sets the current_chunk_group ptr to point to the named chunk group.
    // (nullptr if it doesn't exist). This seems oddly stateful for my
    // tastes but I did it anyway to avoid causing a dangling pointer.
    // (Subsequent world_cache.get_entry may overwrite ptr).
    //
    // I may re-think the design of this program at some point. It
    // used to be quite clean but the quality has eroded with time as
    // I stapled-on more optimizations (and who knows what will happen
    // if I actually switch to Vulkan, my vkFu is weak).
    void set_current_chunk_group(glm::ivec3 group_coord)
    {
        auto& entry = world_cache.get_entry(group_coord);
        current_chunk_group = entry.chunk_group_ptr.get();
    }

  protected:
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
    bool cull_group(glm::ivec3 group_coord, float* squared_dist=nullptr) const
    {
        auto& tr = transforms;
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

    // Return a lambda form of cull_group (passed to the
    // MeshStore/RaycastStore to avoid spending time loading groups
    // that are outside the view frustum).
    auto get_cull_acceptor() const
    {
        auto& self = *this;
        return [&self] (glm::ivec3 coord) { return !self.cull_group(coord); };
    }

    static constexpr int cull = 1, draw_mesh = 2, draw_raycast = 3;
    // Distance culling step. Given the group coordinate and AABB of a
    // chunk, decide which of the above things we should do: cull it,
    // draw it with the mesh algorithm (near the eye), or draw it with
    // the raycast algorithm (far away).
    //
    // This function may be duplicated on the GPU, hence the
    // floor(eye) -- this prevents subtle disagreements due to
    // different FP representations. Also, think twice before changing
    // this function's implementation.
    int decide_chunk(glm::ivec3 group_coord,
                     PackedAABB aabb) const
    {
        auto& tr = transforms;
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

  private:
    // 3D cache of mesh/raycast entries (one entry per chunk group).
    // Really all this is just to dispatch the real work to the
    // abstract functions: see AsyncCache for actual implementation details.
    template <typename StagingT, typename EntryT>
    struct StoreT : public AsyncCache<StagingT, EntryT>
    {
        // Back pointer to RendererLogic.
        RendererLogic* owner;

      public:
        StoreT(RendererLogic* owner_, const AsyncCacheArgs& args) :
            AsyncCache<StagingT, EntryT>(args),
            owner(owner_)
        {

        }

        // Read from disk/memory data for the given chunk group and
        // pass it to the function filling a staging buffer.
        // We promised to fix the dirty flag.
        bool stage(StagingT* staging, glm::ivec3 group_coord) override
        {
            thread_local std::unique_ptr<ViewWorldCache> world_cache_ptr;
            if (world_cache_ptr == nullptr) {
                world_cache_ptr.reset(new ViewWorldCache(owner->world_handle));
            }

            auto& entry = world_cache_ptr->get_entry(group_coord);
            const BinChunkGroup* group_ptr = entry.chunk_group_ptr.get();
            if (group_ptr == nullptr) return false;

            owner->worker_stage(staging, group_ptr);
            return true;
        }

        // Swap from staging buffer to the cache (if ready).
        bool swap_in(StagingT* staging,
                     std::unique_ptr<EntryT>* p_uptr_entry) override
        {
            return owner->swap_in(staging, p_uptr_entry);
        }
    };

    class MeshStore : public StoreT<MeshStaging, MeshEntry>
    {
      public:
        MeshStore(RendererLogic* owner_, const AsyncCacheArgs& args) :
        StoreT<MeshStaging, MeshEntry>(owner_, args) { }
    };

    class RaycastStore : public StoreT<RaycastStaging, RaycastEntry>
    {
      public:
        RaycastStore(RendererLogic* owner_, const AsyncCacheArgs& args) :
        StoreT<RaycastStaging, RaycastEntry>(owner_, args) { }
    };

    // Collect a list of chunk groups (MeshEntry) that are within the
    // view frustum, and close enough for mesh rendering, and send
    // them to the mesh rendering implementation.
    void render_world_mesh_step()
    {
        if (mesh_store == nullptr) {
            static constexpr AsyncCacheArgs args = {
                3,    // modulus
                8,    // associativity
                32,   // staging buffers
                1,    // worker threads
                10,   // condvar timeout in milliseconds
                2,    // Frames elapsed since entry use before eviction allowed
            };
            mesh_store.reset(new MeshStore(this, args));
        }
        std::vector<std::pair<MeshEntry*, glm::ivec3>> entries;
        MeshStore& store = *mesh_store;
        store.begin_frame();

        float squared_thresh = mesh_load_threshold * mesh_load_threshold;

        glm::ivec3 group_coord_low, group_coord_high;
        glm::dvec3 disp(mesh_load_threshold);
        split_coordinate(transforms.eye - disp, &group_coord_low);
        split_coordinate(transforms.eye + disp, &group_coord_high);

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
        store.stage_from_queue(8, get_cull_acceptor());
        store.swap_in_from_staging(8);
    }

    // Collect a list of chunk groups (RaycastEntry) that are within
    // the view frustum and send them to the raycast renderer
    // implementation.
    void render_world_raycast_step()
    {
        if (raycast_store == nullptr) {
            static constexpr AsyncCacheArgs args = {
                4,   // modulus
                12,  // associativity (CHECK THIS IF THERE's FLICKERING)
                128, // staging buffers
                2,   // worker threads
                10,  // condvar timeout in milliseconds
                2,   // Frames elapsed since entry use before eviction allowed
            };
            raycast_store.reset(new RaycastStore(this, args));
        }
        RaycastStore& store = *raycast_store;
        store.begin_frame();

        // Resize the cache to a size suitable for the given render distance.
        auto far_plane = transforms.far_plane;
        uint32_t modulus = uint32_t(std::max(4.0, ceil(far_plane / 128.0)));
        store.set_modulus(modulus, [&] { this->main_thread_wait_idle(); } );

        // Count and limit the number of new chunk groups added to the
        // RaycastStore this frame.
        int remaining_new_chunk_groups_allowed =
            transforms.max_frame_new_chunk_groups;

        // Collect all chunk groups needing to be raycast, and sort from
        // nearest to furthest (reverse painters). This fascilitates
        // the early depth test optimization.
        float squared_thresh = transforms.far_plane;
        squared_thresh *= squared_thresh;
        std::vector<std::pair<float, glm::ivec3>> group_coord_by_depth;

        unsigned drawn_group_count = 0;
        glm::ivec3 group_coord_low, group_coord_high;
        glm::dvec3 disp(transforms.far_plane);
        split_coordinate(transforms.eye - disp, &group_coord_low);
        split_coordinate(transforms.eye + disp, &group_coord_high);

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

        // if (!disable_zcull_sort) {
        if (true) {
            auto lt_depth = [] (const auto& left, const auto& right)
            {
                return left.first < right.first;
            };
            std::sort(group_coord_by_depth.begin(),
                      group_coord_by_depth.end(),
                      lt_depth);
        }

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
        raycast_store->stage_from_queue(80, get_cull_acceptor());
    }

  protected:
    void draw_frame() override final
    {
        transforms = p_camera->get_transforms_vk(); // Use glClipControl.
        begin_frame();
        render_world_mesh_step();
        render_world_raycast_step();
        end_frame();
    }

    void destroy_stores() override final
    {
        mesh_store->stop_wait_threads();
        raycast_store->stop_wait_threads();
        wait_idle();
        mesh_store.reset(nullptr);
        raycast_store.reset(nullptr);
    }
};

} // end namespace myricube

#endif /* !MYRICUBE_RENDERERLOGIC_HH_ */
