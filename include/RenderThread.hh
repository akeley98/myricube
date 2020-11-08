// Handle for a rendering thread (runs OpenGL or whatever renderer in
// a loop on its own thread). Nothing really is happening here except
// storing a handle to a thread and a flag for ordering the thread to
// exit; the real work is in RendererLogic.

#ifndef MYRICUBE_RENDERTHREAD_HH_
#define MYRICUBE_RENDERTHREAD_HH_

#include "myricube.hh"

#include <atomic>
#include <cassert>
#include <memory>
#include <thread>
#include <utility>

#include "window.hh"

namespace myricube {

class RendererLogic;

void render_loop(RendererLogic&);

// Passed to the renderer thread.
struct RenderArgs
{
    // Window to draw to and camera that other threads can modify to
    // control the renderer's camer position.
    // Renderer thread takes shared ownership of these.
    std::shared_ptr<Window> p_window;
    std::shared_ptr<SyncCamera> p_camera;

    // What world to draw.
    WorldHandle world_handle;
};

class RenderThread
{
    // Set to true when the render thread should exit its main loop.
    std::atomic<bool> thread_exit_flag{false};

    // Atomically set by the renderer thread to report stats.
    std::atomic<double> atomic_fps{0};
    std::atomic<double> atomic_frame_time{0};

    // The renderer thread.
    std::thread render_thread;

    friend class RendererLogic;

  public:
    // Start the renderer thread: this thread runs the given
    // RendererLogic factory function and runs the resulting
    // RendererLogic's renderer_loop until ordered to stop. The given
    // RenderArgs is passed through.
    typedef std::shared_ptr<RendererLogic> Factory(RenderThread*, RenderArgs);
    RenderThread(Factory factory, RenderArgs args)
    {
        assert(args.p_window != nullptr);
        assert(args.p_camera != nullptr);

        auto thread_loop = [self=this, factory=factory, args=std::move(args)]
        {
            std::shared_ptr<RendererLogic> ptr = factory(self, std::move(args));
            render_loop(*ptr);
        };
        render_thread = std::thread(std::move(thread_loop));
    }

    // Since the render thread has a pointer back to us, we must not
    // move in memory.
    RenderThread(RenderThread&&) = delete;

    // Implicitly stop the render thread upon destruction.
    ~RenderThread()
    {
        stop();
        render_thread.join();
    }

    // Order the render thread to eventually stop. May safely be
    // called multiple times (i.e. is idempotent).
    void stop()
    {
        thread_exit_flag.store(true);
    }

    // Get the latest fps reported by the renderer thread.
    double get_fps()
    {
        return atomic_fps.get();
    }

    // Get the frame time (in seconds) reported by the renderer thread.
    double get_frame_time()
    {
        return atomic_frame_time.get();
    }
};

} // end namespace myricube

#endif /* !MYRICUBE_RENDERTHREAD_HH_ */
