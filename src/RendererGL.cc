// OpenGL implementation of RendererLogic.
#include "RendererLogic.hh"

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

class RendererGL :
    RendererLogic<MeshEntry, MeshStaging, RaycastEntry, RaycastStaging>
{
    void begin_frame() override
    {

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
};

} // end namespace myricube
