// Functions for drawing the voxel world. As typical for OpenGL, these
// are NOT threadsafe, even for drawing different worlds in parallel.

#ifndef MYRICUBE_RENDERER_HH_
#define MYRICUBE_RENDERER_HH_

namespace myricube {

#include <memory>

class Renderer;
class SyncCamera;

// Not threadsafe for now; todo.
Renderer* new_renderer(const WorldHandle&, std::shared_ptr<SyncCamera>);
void delete_renderer(Renderer*);
void draw_frame(Renderer*);

// Experimental: Trying out 32-bit float depth buffer. This requires
// (in practice) rendering to an offscreen framebuffer.

// Bind the global off-screen framebuffer. Size is given in
// case the framebuffer must be resized or created. We also pass a
// target number of fragments -- if nonzero, the intermediate
// framebuffer may be an (integer) multiple smaller than the actual
// screen size if rendering at full-resolution would exceed the target
// fragment count.
void bind_global_f32_depth_framebuffer(
    int frame_x, int frame_y, int target_fragments=0);

// Blit the contents of the said framebuffer to the window
// framebuffer, and bind the window framebuffer again.
void finish_global_f32_depth_framebuffer(int frame_x, int frame_y);

} // end namespace
#endif /* !MYRICUBE_RENDERER_HH_ */
