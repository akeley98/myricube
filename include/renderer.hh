// Renderer class that spawns a secondary thread running OpenGL code
// for drawing a voxel world.

#ifndef MYRICUBE_RENDERER_HH_
#define MYRICUBE_RENDERER_HH_

namespace myricube {

#include <memory>

class Renderer;
class SyncCamera;
class Window;

// TODO: Multiple Renderers won't work correctly due to global state.

// Start up a thread drawing the given voxel world. The OpenGL context
// of the given Window is used, and the camera can be controlled from
// other threads through the shared SyncCamera.
Renderer* new_renderer(
    std::shared_ptr<Window>,
    WorldHandle,
    std::shared_ptr<SyncCamera>);

// Stop the renderer thread, then clean up its resources.
void delete_renderer(Renderer*);

// Get some RAII going for Renderer.
struct RendererDeleter
{
    void operator() (Renderer* victim) { delete_renderer(victim); }
};
using UPtrRenderer = std::unique_ptr<Renderer, RendererDeleter>;

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
