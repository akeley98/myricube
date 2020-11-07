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
// functions". Also, when I say "lists of chunk group", what I really
// mean is lists of the MeshEntry/RaycastEntry template parameters,
// which should hold the 3D-textures or whatever needed to actually
// render a chunk group.

#ifndef MYRICUBE_RENDERERLOGIC_HH_
#define MYRICUBE_RENDERERLOGIC_HH_

#include "myricube.hh"

#include <cassert>
#include <thread>
#include <vector>

#include "AsyncCache.hh"
#include "camera.hh"
#include "RenderThread.hh"
#include "voxels.hh"

// Need AABB and Mesh header files aabb.hh and mesh.hh.

namespace myricube {

// MeshEntry and MeshShading are the EntryT and StagingT type
// parameters for the AsyncCache used for mesh rendering
// (i.e. MeshEntry should be the API data needed to actually render a
// chunk group with the implemented mesh renderer). RaycastEntry and
// RaycastStaging are EntryT and StagingT for raycast rendering.
//
// Throughout, "main thread" refers to the thread that created this
// RendererLogic and is running its render_loop.
template <MeshEntry, MeshStaging, RaycastEntry, RaycastStaging>
class RendererLogic
{
    // Back pointer to the RenderThread that instantiated this
    // RendererLogic.
    RenderThread* const p_back = nullptr;
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
        p_back(p_back_),
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
    // draw them using the mesh rendering strategy.
    virtual void
    draw_mesh_entries(
        const std::vector<std::pair<MeshEntry*, glm::ivec3>>&) = 0;

    // Called by the main thread after draw_mesh_entries.
    //
    // Given a vector of chunk groups (represented as pairs of their
    // RaycastEntry and group coordinate), draw them using the raycast
    // rendering strategy.
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
    // Return true if successful, false if we need to retry later.
    virtual bool swap_in(MeshStaging*, std::unique_ptr<MeshEntry>*) = 0;

    // Same, but for raycasting.
    virtual bool swap_in(RaycastStaging*, std::unique_ptr<RaycastEntry>*) = 0;

    // NOTE: I really hope the above stays up-to-date if I change the
    // AsyncCache...
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
    bool cull_group(glm::ivec3 group_coord, float* squared_dist=nullptr)
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
};

} // end namespace myricube

#endif /* !MYRICUBE_RENDERERLOGIC_HH_ */
