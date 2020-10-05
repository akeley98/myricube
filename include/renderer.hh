// Renderer class that spawns a secondary thread running OpenGL code
// for drawing a voxel world.

#ifndef MYRICUBE_RENDERER_HH_
#define MYRICUBE_RENDERER_HH_

namespace myricube {

#include <memory>

class Renderer;
class SyncCamera;
class Window;
class WorldHandle;

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

// Get (approximate) FPS and frame time reported by the Renderer.
double get_fps(const Renderer&);
double get_frame_time(const Renderer&);

// Get said frame time in integer milliseconds.
inline int get_frame_time_ms(const Renderer& renderer)
{
    return int(get_frame_time(renderer) * 1000);
}

} // end namespace
#endif /* !MYRICUBE_RENDERER_HH_ */
