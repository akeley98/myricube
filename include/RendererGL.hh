#ifndef MYRICUBE_RENDERERGL_HH_
#define MYRICUBE_RENDERERGL_HH_

#include "myricube.hh"
#include "RenderThread.hh"

#include <memory>

namespace myricube {

std::shared_ptr<RendererBase> RendererGL_Factory(RenderThread*, RenderArgs);

// Alternate version that does not use compute shaders.
std::shared_ptr<RendererBase> RendererGL_Factory_glTex(
    RenderThread*, RenderArgs);

} // end namespace myricube

#endif /* !MYRICUBE_RENDERERGL_HH_ */

