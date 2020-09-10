// Functions for drawing the voxel world. As typical for OpenGL, these
// are NOT threadsafe, even for drawing different worlds in parallel.

#ifndef MYRICUBE_RENDERER_HH_
#define MYRICUBE_RENDERER_HH_

namespace myricube {

class RaycastStore
{
  public:
    virtual ~RaycastStore() = default;
};
class MeshStore;

RaycastStore* new_raycast_store();
MeshStore* new_mesh_store();
void delete_raycast_store(RaycastStore*);
void delete_mesh_store(MeshStore*);

class Camera;
class VoxelWorld;

// Render, to the current framebuffer, chunks near the camera using
// the conventional mesh-based algorithm.
void render_world_mesh_step(VoxelWorld&, Camera&);

// Render, to the current framebuffer, chunks around the camera
// using the AABB-raycast algorithm. Chunks that are near
// enough to have been drawn using the mesh algorithm will
// not be re-drawn.
void render_world_raycast_step(VoxelWorld&, Camera&);

// Render the background. ("skybox")
void render_background(Camera&);

void gl_clear();

// Wrapper for glViewport
void viewport(int x, int y);

// glEnables and stuff.
void gl_first_time_setup();

// Experimental: Trying out 32-bit float depth buffer. This requires
// (in practice) rendering to an offscreen framebuffer.

// Bind the global off-screen framebuffer. Screen size is given in
// case the framebuffer must be resized or created.
void bind_global_f32_depth_framebuffer(int screen_x, int screen_y);

// Blit the contents of the said framebuffer to the window
// framebuffer, and bind the window framebuffer again.
void finish_global_f32_depth_framebuffer(int screen_x, int screen_y);

} // end namespace
#endif /* !MYRICUBE_RENDERER_HH_ */
