// Dummy implementation of RendererGL
#define GLFW_INCLUDE_NONE
#include "myricube.hh"
#include "RendererLogic.hh"

namespace myricube {

std::shared_ptr<RendererBase> RendererGL_Factory(
    RenderThread*,
    RenderArgs)
{
    panic("OpenGL renderer not compiled in (see src/DummyGL.cc)");
    return nullptr;
}

}
