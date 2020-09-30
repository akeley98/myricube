// Functions for drawing the voxel world. As typical for OpenGL, these
// are NOT threadsafe, even for drawing different worlds in parallel.

#ifndef MYRICUBE_RENDERER_HH_
#define MYRICUBE_RENDERER_HH_

namespace myricube {

class Renderer;

class RaycastStore;
class MeshStore;
class Camera;
class WorldHandle;

// Not threadsafe for now; todo.
Renderer* new_renderer(const WorldHandle&, Camera*);
void delete_renderer(Renderer*);
void draw_frame(Renderer*);

RaycastStore* new_raycast_store();
MeshStore* new_mesh_store();
void delete_raycast_store(RaycastStore*);
void delete_mesh_store(MeshStore*);

// Experimental: Trying out 32-bit float depth buffer. This requires
// (in practice) rendering to an offscreen framebuffer.

// Bind the global off-screen framebuffer. Screen size is given in
// case the framebuffer must be resized or created. We also pass a
// target number of fragments -- if nonzero, the intermediate
// framebuffer may be an (integer) multiple smaller than the actual
// screen size if rendering at full-resolution would exceed the target
// fragment count.
void bind_global_f32_depth_framebuffer(
    int screen_x, int screen_y, int target_fragments=0);

// Blit the contents of the said framebuffer to the window
// framebuffer, and bind the window framebuffer again.
void finish_global_f32_depth_framebuffer(int screen_x, int screen_y);

} // end namespace
#endif /* !MYRICUBE_RENDERER_HH_ */
