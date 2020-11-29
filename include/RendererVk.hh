#ifndef MYRICUBE_RENDERERVK_HH_
#define MYRICUBE_RENDERERVK_HH_

#include "myricube.hh"
#include "RenderThread.hh"

#include <memory>

namespace myricube {

extern bool vulkan_compiled_in;
std::shared_ptr<RendererBase> RendererVk_Factory(RenderThread*, RenderArgs);

} // end namespace myricube

#endif /* !MYRICUBE_RENDERERGL_HH_ */

