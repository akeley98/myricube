// Handle for a rendering thread (runs OpenGL or whatever renderer in
// a loop on its own thread). Nothing really is happening here except
// implementing the top-level render loop (calls draw_frame in a
// loop), and checking for an atomic flag that tells the loop to
// stop. the real work is in RendererLogic.

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

class RenderThread;

// Non-templatized base class for RendererLogic.
class RendererBase
{
    // Back pointer to the RenderThread that instantiated this
    // RendererLogic.
    RenderThread* const p_back = nullptr;
    friend class RenderThread;

  protected:
    RendererBase(RenderThread* p_back_) : p_back(p_back_) { }
    virtual void draw_frame() = 0;
    virtual void destroy_stores() = 0;

  private:
    inline void render_loop();
};

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

    friend class RendererBase;

  public:
    // Start the renderer thread: this thread runs the given
    // RendererLogic factory function and runs the resulting
    // RendererLogic's renderer_loop until ordered to stop. The given
    // RenderArgs is passed through.
    typedef std::shared_ptr<RendererBase> Factory(RenderThread*, RenderArgs);
    RenderThread(Factory factory, RenderArgs args)
    {
        assert(args.p_window != nullptr);
        assert(args.p_camera != nullptr);

        auto thread_loop = [self=this, factory=factory, args=std::move(args)]
        {
            std::shared_ptr<RendererBase> ptr = factory(self, std::move(args));
            ptr->render_loop();
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
    double get_fps() const
    {
        return atomic_fps.load();
    }

    // Get the frame time (in seconds) reported by the renderer thread.
    double get_frame_time() const
    {
        return atomic_frame_time.load();
    }
};

void RendererBase::render_loop()
{
    double previous_update = glfwGetTime();
    double previous_fps_update = glfwGetTime();
    int frames = 0;
    double frame_time = 0;

    // Draw frames in a loop until the renderer is ordered to stop,
    // calculating (CPU-side) FPS and frame time.
    while (!p_back->thread_exit_flag.load()) {
        draw_frame();

        // Calculate (approximate) frame time. Require at least 2 ms
        // between frames (workaround to system freeze bug).
        double now, dt;
        do {
            now = glfwGetTime();
            dt = now - previous_update;
        } while (dt < 0.002);
        previous_update = now;

        // Update FPS and frame time. Periodically report back to the
        // RenderThread.
        ++frames;
        if (dt > frame_time) frame_time = dt;

        if (now - previous_fps_update >= 0.5) {
            p_back->atomic_fps.store(frames / (now - previous_fps_update));
            p_back->atomic_frame_time.store(frame_time);

            previous_fps_update = now;
            frames = 0;
            frame_time = 0;
        }
    }

    // DON'T REMOVE THIS.
    destroy_stores();
}

} // end namespace myricube

#endif /* !MYRICUBE_RENDERTHREAD_HH_ */
