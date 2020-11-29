// Dummy implementation of RendererVk
#define GLFW_INCLUDE_NONE
#include "myricube.hh"
#include "RendererLogic.hh"

namespace myricube {

bool vulkan_compiled_in = false;

std::shared_ptr<RendererBase> RendererVk_Factory(
    RenderThread*,
    RenderArgs)
{
    panic("Vulkan renderer not compiled in (see src/DummyVk.cc)");
    return nullptr;
}

}
