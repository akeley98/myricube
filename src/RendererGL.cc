// OpenGL implementation of RendererLogic.
#include "RendererLogic.hh"

#include "glad/glad.h"

namespace {

struct MeshEntry
{

};

struct MeshStaging
{

};

struct RaycastEntry
{

};

struct RaycastStaging
{

};

} // end anonymous namespace.

namespace myricube {

struct RendererGL :
    RendererLogic<MeshEntry, MeshStaging, RaycastEntry, RaycastStaging>
{
    RendererGL(RenderThread* thread, RenderArgs args) :
        RendererLogic<MeshEntry, MeshStaging, RaycastEntry, RaycastStaging>
        (thread, args)
    {
        p_window->gl_make_current();
        glClearColor(1, 1, 0, 1);
    }

    void begin_frame() override
    {
        glClear(GL_COLOR_BUFFER_BIT);
    }

    void
    draw_mesh_entries(
        const std::vector<std::pair<MeshEntry*, glm::ivec3>>&) override
    {

    }

    void
    draw_raycast_entries(
        const std::vector<std::pair<RaycastEntry*, glm::ivec3>>&) override
    {

    }

    void end_frame() override
    {
        p_window->swap_buffers();
    }

    void worker_stage(MeshStaging*, const BinChunkGroup*) override
    {

    }

    void worker_stage(RaycastStaging*, const BinChunkGroup*) override
    {

    }

    bool swap_in(MeshStaging*, std::unique_ptr<MeshEntry>*) override
    {
        return true;
    }

    bool swap_in(RaycastStaging*, std::unique_ptr<RaycastEntry>*) override
    {
        return true;
    }

    void wait_idle() override
    {
        glFinish();
    }
};

std::shared_ptr<RendererBase> RendererGL_Factory(
    RenderThread* thread,
    RenderArgs args)
{
    return std::shared_ptr<RendererBase>(new RendererGL(thread, args));
}

} // end namespace myricube
