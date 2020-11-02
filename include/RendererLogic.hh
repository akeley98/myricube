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
// the loop.
//
// The main purpose of this class is to manage at a high level the GPU
// caches of chunk groups (see AsyncCache), load chunk groups from
// disk/memory and add them to the caches as-needed, and extract lists
// of chunk groups from the caches (culled based on camera location)
// and pass them on to the subclass's functions to render them.
//
// In the above, when I say "chunk group", what I really mean is the
// MeshEntry/RaycastEntry template parameters, which should hold the
// 3D-textures or whatever needed to actually render a chunk group.

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
    RenderThread* back_pointer;
  protected:
    // Window, camera, world arguments used to launch the render thread.
    std::shared_ptr<Window> p_window;
    std::shared_ptr<SyncCamera> p_camera;
    WorldHandle world_handle;

    // Constructed from world_handle to efficiently view chunk groups
    // of the world. Access only through set_current_chunk_group.
    ViewWorldCache world_cache;

    // A bit hacky, but this is the last chunk group requested from
    // world_cache. Needed because new accesses invalidate
    // previously-returned pointers from world_cache.
    const BinChunkGroup* current_chunk_group = nullptr;

    // Updated from SyncCamera every frame.
    CameraTransforms transforms;

  private:
    AsyncCache<MeshEntry, MeshStaging> mesh_cache;
    AsyncCache<RaycastEntry, RaycastStaging> raycast_cache;

  protected:
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
    // group.
    virtual void worker_stage(MeshStaging*, const BinChunkGroup*) = 0;

    // Like the above, but for raycasting. Same warning applies.
    virtual void worker_stage(RaycastStaging*, const BinChunkGroup*) = 0;

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
};

} // end namespace myricube

#endif /* !MYRICUBE_RENDERERLOGIC_HH_ */
