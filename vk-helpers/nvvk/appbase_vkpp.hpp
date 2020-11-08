/* Copyright (c) 2014-2018, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#pragma once

#include <vulkan/vulkan.hpp>

#include "imgui.h"
#include "imgui_helper.h"
#include "imgui_impl_vk.h"
#include "nvh/camerainertia.hpp"
#include "nvh/cameramanipulator.hpp"
#include "swapchain_vk.hpp"

#ifdef LINUX
#include <unistd.h>
#endif

#ifdef _WIN32
#define GLFW_EXPOSE_NATIVE_WIN32
#endif
#include "GLFW/glfw3.h"
#include "GLFW/glfw3native.h"

#include <cmath>
#include <vector>

namespace nvvk {
static float const s_keyTau    = 0.10f;
static float const s_cameraTau = 0.03f;
static float const s_moveStep  = 0.2f;

//--------------------------------------------------------------------------------------------------
/**

# class nvvk::AppBase

The framework comes with a few `App???` classes, these can serve as base classes for various samples.
They might differ a bit in setup and functionality, but in principle aid the setup of context and window,
as well as some common event processing.

The nvvk::AppBase serves as the base class for many ray tracing examples and makes use of the Vulkan C++ API (`vulkan.hpp`).
It does the basics for Vulkan, by holding a reference to the instance and device, but also comes with optional default setups 
for the render passes and the swapchain.

## Usage

An example will derive from this class:

~~~~ C++
class MyExample : public AppBase 
{
};
~~~~

## Setup

In the `main()` of the example, after creating the Vulkan instance and device, call `setup()`.
This will hold the VkInstance, VkDevice, VkPhysicalDevice and create the VkQueue, VkPool, plus it
will initialize all Vulkan extensions for the C++ API (vulkan.hpp).

Prior to calling setup, if you are using the `nvvk::Context` class to create and initalize Vulkan instances,
you may want to create a VkSurfaceKHR from the window (glfw for example) and call `setGCTQueueWithPresent()`.
This will make sure the Queue indices are adapted.

Creating the surface, will actually create the swapchain for displaying. Arguments are
width and height, color and depth format, and vsync on/off. Defaults will create the best format.

Before creating the framebuffers to display the results, a depth buffer has to be create which will
be used by all framebuffers. Then we can create the framebuffer.

**Note**: the imageView(s) are part of the swapchain. 

There is also a 'default renderpass', which is a color/depth, clear both buffers.

If the application is using Dear ImGui, there are convenient functions for initializing it and
setting the callbacks (glfw). The first one to call is `initGUI(0)`, where the argument is the subpass
where it will be use. Default is 0, but if the application creates a renderpass with multi-sampling and
resolves in the second subpass, this makes it possible.

Then call `setupGlfwCallbacks(window)` to have all the window callback: key, mouse, window resizing.
By default AppBase will handle resizing of the window. It will recreate the images and framebuffers,
but a sample may need to overload that function, or to be aware of that change, therefore can overload
`onResize(width, height)`.

Last setup for ImGui is to call `ImGui_ImplGlfw_InitForVulkan(window, true)`, where true is for the 
callbacks for Imgui.

**Note**: All the methods are virtual and can be overloaded if they are not doing the typical setup. 

~~~~ C++
MyExample example;

const vk::SurfaceKHR surface = example.getVkSurface(vkctx.m_instance, window);
vkctx.setGCTQueueWithPresent(surface);

example.setup(vkctx.m_instance, vkctx.m_device, vkctx.m_physicalDevice, vkctx.m_queueGCT.familyIndex);
example.createSurface(surface, SAMPLE_SIZE_WIDTH, SAMPLE_SIZE_HEIGHT);
example.createDepthBuffer();
example.createFrameBuffers();
example.createRenderPass();
example.initGUI(0);
example.setupGlfwCallbacks(window);

ImGui_ImplGlfw_InitForVulkan(window, true);
~~~~

## Drawing loop

The drawing loop in the main() is the typicall loop you will find in glfw examples. Note that
AppBase has a convenient function to tell if the window is minimize, therefore not doing any 
work and contain a sleep(), so the CPU is not going crazy. 


~~~~ C++
// Window system loop
while(!glfwWindowShouldClose(window))
{
  glfwPollEvents();
  if(example.isMinimized())
    continue;

  example.display();  // infinitely drawing
}
~~~~

## Display

A typical display() function will need the following: 

* Acquiring the next image: `prepareFrame()`
* Get the command buffer for the frame. There are n command buffers equal to the number of in-flight frames.
* Clearing values
* Start rendering pass
* Drawing
* End rendering
* Submitting frame to display

~~~~ C++
void MyExample::display()
{
  // Acquire 
  prepareFrame();

  // Command buffer for current frame
  const vk::CommandBuffer& cmdBuff = m_commandBuffers[getCurFrame()];
  cmdBuff.begin({vk::CommandBufferUsageFlagBits::eOneTimeSubmit});

  // Clearing values
  vk::ClearValue clearValues[2];
  clearValues[0].setColor(std::array<float, 4>({0.1f, 0.1f, 0.4f, 0.f}));
  clearValues[1].setDepthStencil({1.0f, 0});

  // Begin rendering
  vk::RenderPassBeginInfo renderPassBeginInfo{m_renderPass, m_framebuffers[getCurFrame()], {{}, m_size}, 2, clearValues};
  cmdBuff.beginRenderPass(renderPassBeginInfo, vk::SubpassContents::eInline);
  
  // .. draw scene ...

  // Draw UI
  ImGui::RenderDrawDataVK(cmdBuff, ImGui::GetDrawData());

  // End rendering
  cmdBuff.endRenderPass();

  // End of the frame and present the one which is ready
  cmdBuff.end();
  submitFrame();
}
~~~~~

## Closing

Finally, all resources can be destroyed by calling `destroy()` at the end of main().

~~~~ C++
  example.destroy();
~~~~

*/

class AppBase
{
public:
  AppBase()          = default;
  virtual ~AppBase() = default;

  virtual void onResize(int /*w*/, int /*h*/){};  // To implement when the size of the window change

  //--------------------------------------------------------------------------------------------------
  // Setup the low level Vulkan for various operations
  //
  virtual void setup(const vk::Instance& instance, const vk::Device& device, const vk::PhysicalDevice& physicalDevice, uint32_t graphicsQueueIndex)
  {
    // Initialize function pointers
    vk::DynamicLoader         dl;
    PFN_vkGetInstanceProcAddr vkGetInstanceProcAddr =
        dl.getProcAddress<PFN_vkGetInstanceProcAddr>("vkGetInstanceProcAddr");
    VULKAN_HPP_DEFAULT_DISPATCHER.init(vkGetInstanceProcAddr);
    VULKAN_HPP_DEFAULT_DISPATCHER.init(instance);
    VULKAN_HPP_DEFAULT_DISPATCHER.init(device);

    m_instance           = instance;
    m_device             = device;
    m_physicalDevice     = physicalDevice;
    m_graphicsQueueIndex = graphicsQueueIndex;
    m_queue              = m_device.getQueue(m_graphicsQueueIndex, 0);
    m_cmdPool = m_device.createCommandPool({vk::CommandPoolCreateFlagBits::eResetCommandBuffer, graphicsQueueIndex});
    m_pipelineCache = device.createPipelineCache(vk::PipelineCacheCreateInfo());
  }

  //--------------------------------------------------------------------------------------------------
  // To call on exit
  //
  virtual void destroy()
  {
    m_device.waitIdle();

    if(ImGui::GetCurrentContext() != nullptr)
    {
      ImGui::ShutdownVK();
      ImGui::DestroyContext();
    }

    m_device.destroyRenderPass(m_renderPass);
    m_device.destroyImageView(m_depthView);
    m_device.destroyImage(m_depthImage);
    m_device.freeMemory(m_depthMemory);
    m_device.destroyPipelineCache(m_pipelineCache);

    for(uint32_t i = 0; i < m_swapChain.getImageCount(); i++)
    {
      m_device.destroyFence(m_waitFences[i]);
      m_device.destroyFramebuffer(m_framebuffers[i]);
      m_device.freeCommandBuffers(m_cmdPool, m_commandBuffers[i]);
    }
    m_swapChain.deinit();

    m_device.destroyCommandPool(m_cmdPool);
    if(m_surface)
      m_instance.destroySurfaceKHR(m_surface);
  }


  //--------------------------------------------------------------------------------------------------
  // Return the surface "screen" for the display
  //
  VkSurfaceKHR getVkSurface(const vk::Instance& instance, GLFWwindow* window)
  {
    assert(instance);
    m_window = window;

    VkSurfaceKHR surface{};
    VkResult     err = glfwCreateWindowSurface(instance, window, nullptr, &surface);
    if(err != VK_SUCCESS)
    {
      assert(!"Failed to create a Window surface");
    }
    m_surface = surface;

    return surface;
  }


  //--------------------------------------------------------------------------------------------------
  // Creating the surface for rendering
  //
  void createSurface(const vk::SurfaceKHR& surface,
                     uint32_t              width,
                     uint32_t              height,
                     vk::Format            colorFormat = vk::Format::eB8G8R8A8Unorm,
                     vk::Format            depthFormat = vk::Format::eD32SfloatS8Uint,
                     bool                  vsync       = false)
  {
    m_size        = vk::Extent2D(width, height);
    m_depthFormat = depthFormat;
    m_colorFormat = colorFormat;
    m_vsync       = vsync;

    m_swapChain.init(m_device, m_physicalDevice, m_queue, m_graphicsQueueIndex, surface, static_cast<VkFormat>(colorFormat));
    m_swapChain.update(m_size.width, m_size.height, vsync);
    m_colorFormat = static_cast<vk::Format>(m_swapChain.getFormat());

    // Create Synchronization Primitives
    m_waitFences.resize(m_swapChain.getImageCount());
    for(auto& fence : m_waitFences)
    {
      fence = m_device.createFence({vk::FenceCreateFlagBits::eSignaled});
    }

    // Command buffers store a reference to the frame buffer inside their render pass info
    // so for static usage without having to rebuild them each frame, we use one per frame buffer
    m_commandBuffers =
        m_device.allocateCommandBuffers({m_cmdPool, vk::CommandBufferLevel::ePrimary, m_swapChain.getImageCount()});

#ifdef _DEBUG
    for(size_t i = 0; i < m_commandBuffers.size(); i++)
    {
      std::string name = std::string("AppBase") + std::to_string(i);
      m_device.setDebugUtilsObjectNameEXT(
          {vk::ObjectType::eCommandBuffer, reinterpret_cast<const uint64_t&>(m_commandBuffers[i]), name.c_str()});
    }
#endif  // _DEBUG

    // Setup camera
    CameraManip.setWindowSize(m_size.width, m_size.height);
  }

  //--------------------------------------------------------------------------------------------------
  // Create the framebuffers in which the image will be rendered
  // - Swapchain need to be created before calling this
  //
  virtual void createFrameBuffers()
  {
    // Recreate the frame buffers
    for(auto framebuffer : m_framebuffers)
    {
      m_device.destroyFramebuffer(framebuffer);
    }

    // Array of attachment (color, depth)
    std::array<vk::ImageView, 2> attachments;

    // Create frame buffers for every swap chain image
    vk::FramebufferCreateInfo framebufferCreateInfo;
    framebufferCreateInfo.renderPass      = m_renderPass;
    framebufferCreateInfo.attachmentCount = 2;
    framebufferCreateInfo.width           = m_size.width;
    framebufferCreateInfo.height          = m_size.height;
    framebufferCreateInfo.layers          = 1;
    framebufferCreateInfo.pAttachments    = attachments.data();

    // Create frame buffers for every swap chain image
    m_framebuffers.resize(m_swapChain.getImageCount());
    for(uint32_t i = 0; i < m_swapChain.getImageCount(); i++)
    {
      attachments[0]    = m_swapChain.getImageView(i);
      attachments[1]    = m_depthView;
      m_framebuffers[i] = m_device.createFramebuffer(framebufferCreateInfo);
    }


#ifdef _DEBUG
    for(size_t i = 0; i < m_framebuffers.size(); i++)
    {
      std::string name = std::string("AppBase") + std::to_string(i);
      m_device.setDebugUtilsObjectNameEXT(
          {vk::ObjectType::eFramebuffer, reinterpret_cast<const uint64_t&>(m_framebuffers[i]), name.c_str()});
    }
#endif  // _DEBUG
  }

  //--------------------------------------------------------------------------------------------------
  // Creating a default render pass, very simple one.
  // Other examples will mostly override this one.
  //
  virtual void createRenderPass()
  {
    if(m_renderPass)
    {
      m_device.destroyRenderPass(m_renderPass);
    }

    std::array<vk::AttachmentDescription, 2> attachments;
    // Color attachment
    attachments[0].setFormat(m_colorFormat);
    attachments[0].setLoadOp(vk::AttachmentLoadOp::eClear);
    attachments[0].setFinalLayout(vk::ImageLayout::ePresentSrcKHR);

    // Depth attachment
    attachments[1].setFormat(m_depthFormat);
    attachments[1].setLoadOp(vk::AttachmentLoadOp::eClear);
    attachments[1].setStencilLoadOp(vk::AttachmentLoadOp::eClear);
    attachments[1].setFinalLayout(vk::ImageLayout::eDepthStencilAttachmentOptimal);

    // One color, one depth
    const vk::AttachmentReference colorReference{0, vk::ImageLayout::eColorAttachmentOptimal};
    const vk::AttachmentReference depthReference{1, vk::ImageLayout::eDepthStencilAttachmentOptimal};

    std::array<vk::SubpassDependency, 1> subpassDependencies;
    // Transition from final to initial (VK_SUBPASS_EXTERNAL refers to all commands executed outside of the actual renderpass)
    subpassDependencies[0].setSrcSubpass(VK_SUBPASS_EXTERNAL);
    subpassDependencies[0].setDstSubpass(0);
    subpassDependencies[0].setSrcStageMask(vk::PipelineStageFlagBits::eBottomOfPipe);
    subpassDependencies[0].setDstStageMask(vk::PipelineStageFlagBits::eColorAttachmentOutput);
    subpassDependencies[0].setSrcAccessMask(vk::AccessFlagBits::eMemoryRead);
    subpassDependencies[0].setDstAccessMask(vk::AccessFlagBits::eColorAttachmentRead | vk::AccessFlagBits::eColorAttachmentWrite);
    subpassDependencies[0].setDependencyFlags(vk::DependencyFlagBits::eByRegion);

    vk::SubpassDescription subpassDescription;
    subpassDescription.setPipelineBindPoint(vk::PipelineBindPoint::eGraphics);
    subpassDescription.setColorAttachmentCount(1);
    subpassDescription.setPColorAttachments(&colorReference);
    subpassDescription.setPDepthStencilAttachment(&depthReference);

    vk::RenderPassCreateInfo renderPassInfo;
    renderPassInfo.setAttachmentCount(static_cast<uint32_t>(attachments.size()));
    renderPassInfo.setPAttachments(attachments.data());
    renderPassInfo.setSubpassCount(1);
    renderPassInfo.setPSubpasses(&subpassDescription);
    renderPassInfo.setDependencyCount(static_cast<uint32_t>(subpassDependencies.size()));
    renderPassInfo.setPDependencies(subpassDependencies.data());

    m_renderPass = m_device.createRenderPass(renderPassInfo);

#ifdef _DEBUG
    m_device.setDebugUtilsObjectNameEXT(
        {vk::ObjectType::eRenderPass, reinterpret_cast<const uint64_t&>(m_renderPass), "AppBase"});
#endif  // _DEBUG
  }

  //--------------------------------------------------------------------------------------------------
  // Creating an image to be used as depth buffer
  //
  virtual void createDepthBuffer()
  {
    if(m_depthView)
      m_device.destroyImageView(m_depthView);
    if(m_depthImage)
      m_device.destroyImage(m_depthImage);
    if(m_depthMemory)
      m_device.freeMemory(m_depthMemory);

    // Depth information
    const vk::ImageAspectFlags aspect = vk::ImageAspectFlagBits::eDepth | vk::ImageAspectFlagBits::eStencil;
    vk::ImageCreateInfo        depthStencilCreateInfo;
    depthStencilCreateInfo.setImageType(vk::ImageType::e2D);
    depthStencilCreateInfo.setExtent(vk::Extent3D{m_size.width, m_size.height, 1});
    depthStencilCreateInfo.setFormat(m_depthFormat);
    depthStencilCreateInfo.setMipLevels(1);
    depthStencilCreateInfo.setArrayLayers(1);
    depthStencilCreateInfo.setUsage(vk::ImageUsageFlagBits::eDepthStencilAttachment | vk::ImageUsageFlagBits::eTransferSrc);
    // Create the depth image
    m_depthImage = m_device.createImage(depthStencilCreateInfo);

    // Allocate the memory
    const vk::MemoryRequirements memReqs = m_device.getImageMemoryRequirements(m_depthImage);
    vk::MemoryAllocateInfo       memAllocInfo;
    memAllocInfo.allocationSize  = memReqs.size;
    memAllocInfo.memoryTypeIndex = getMemoryType(memReqs.memoryTypeBits, vk::MemoryPropertyFlagBits::eDeviceLocal);
    m_depthMemory                = m_device.allocateMemory(memAllocInfo);

    // Bind image and memory
    m_device.bindImageMemory(m_depthImage, m_depthMemory, 0);

    // Create an image barrier to change the layout from undefined to DepthStencilAttachmentOptimal
    vk::CommandBuffer             cmdBuffer;
    vk::CommandBufferAllocateInfo cmdBufAllocateInfo;
    cmdBufAllocateInfo.commandPool        = m_cmdPool;
    cmdBufAllocateInfo.level              = vk::CommandBufferLevel::ePrimary;
    cmdBufAllocateInfo.commandBufferCount = 1;
    cmdBuffer                             = m_device.allocateCommandBuffers(cmdBufAllocateInfo)[0];
    cmdBuffer.begin(vk::CommandBufferBeginInfo{vk::CommandBufferUsageFlagBits::eOneTimeSubmit});

    // Put barrier on top, Put barrier inside setup command buffer
    vk::ImageSubresourceRange subresourceRange;
    subresourceRange.aspectMask = aspect;
    subresourceRange.levelCount = 1;
    subresourceRange.layerCount = 1;
    vk::ImageMemoryBarrier imageMemoryBarrier;
    imageMemoryBarrier.oldLayout               = vk::ImageLayout::eUndefined;
    imageMemoryBarrier.newLayout               = vk::ImageLayout::eDepthStencilAttachmentOptimal;
    imageMemoryBarrier.image                   = m_depthImage;
    imageMemoryBarrier.subresourceRange        = subresourceRange;
    imageMemoryBarrier.srcAccessMask           = vk::AccessFlags();
    imageMemoryBarrier.dstAccessMask           = vk::AccessFlagBits::eDepthStencilAttachmentWrite;
    const vk::PipelineStageFlags srcStageMask  = vk::PipelineStageFlagBits::eTopOfPipe;
    const vk::PipelineStageFlags destStageMask = vk::PipelineStageFlagBits::eEarlyFragmentTests;

    cmdBuffer.pipelineBarrier(srcStageMask, destStageMask, vk::DependencyFlags(), nullptr, nullptr, imageMemoryBarrier);
    cmdBuffer.end();
    m_queue.submit(vk::SubmitInfo{0, nullptr, nullptr, 1, &cmdBuffer}, vk::Fence());
    m_queue.waitIdle();
    m_device.freeCommandBuffers(m_cmdPool, cmdBuffer);

    // Setting up the view
    vk::ImageViewCreateInfo depthStencilView;
    depthStencilView.setViewType(vk::ImageViewType::e2D);
    depthStencilView.setFormat(m_depthFormat);
    depthStencilView.setSubresourceRange({aspect, 0, 1, 0, 1});
    depthStencilView.setImage(m_depthImage);
    m_depthView = m_device.createImageView(depthStencilView);
  }

  //--------------------------------------------------------------------------------------------------
  // Convenient function to call before rendering
  //
  void prepareFrame()
  {
    // Acquire the next image from the swap chain
    if(!m_swapChain.acquire())
    {
      assert(!"This shouldn't happen");
    }

    // Use a fence to wait until the command buffer has finished execution before using it again
    uint32_t imageIndex = m_swapChain.getActiveImageIndex();
    while(m_device.waitForFences(m_waitFences[imageIndex], VK_TRUE, 10000) == vk::Result::eTimeout)
    {
    }
  }

  //--------------------------------------------------------------------------------------------------
  // Convenient function to call for submitting the rendering command
  //
  void submitFrame()
  {
    uint32_t imageIndex = m_swapChain.getActiveImageIndex();
    m_device.resetFences(m_waitFences[imageIndex]);

    // In case of using NVLINK
    const uint32_t                deviceMask  = m_useNvlink ? 0b0000'0011 : 0b0000'0001;
    const std::array<uint32_t, 2> deviceIndex = {0, 1};

    vk::DeviceGroupSubmitInfo deviceGroupSubmitInfo;
    deviceGroupSubmitInfo.setWaitSemaphoreCount(1);
    deviceGroupSubmitInfo.setCommandBufferCount(1);
    deviceGroupSubmitInfo.setPCommandBufferDeviceMasks(&deviceMask);
    deviceGroupSubmitInfo.setSignalSemaphoreCount(m_useNvlink ? 2 : 1);
    deviceGroupSubmitInfo.setPSignalSemaphoreDeviceIndices(deviceIndex.data());
    deviceGroupSubmitInfo.setPWaitSemaphoreDeviceIndices(deviceIndex.data());

    vk::Semaphore semaphoreRead  = m_swapChain.getActiveReadSemaphore();
    vk::Semaphore semaphoreWrite = m_swapChain.getActiveWrittenSemaphore();

    // Pipeline stage at which the queue submission will wait (via pWaitSemaphores)
    const vk::PipelineStageFlags waitStageMask = vk::PipelineStageFlagBits::eColorAttachmentOutput;
    // The submit info structure specifies a command buffer queue submission batch
    vk::SubmitInfo submitInfo;
    submitInfo.setPWaitDstStageMask(&waitStageMask);  // Pointer to the list of pipeline stages that the semaphore waits will occur at
    submitInfo.setPWaitSemaphores(&semaphoreRead);  // Semaphore(s) to wait upon before the submitted command buffer starts executing
    submitInfo.setWaitSemaphoreCount(1);               // One wait semaphore
    submitInfo.setPSignalSemaphores(&semaphoreWrite);  // Semaphore(s) to be signaled when command buffers have completed
    submitInfo.setSignalSemaphoreCount(1);             // One signal semaphore
    submitInfo.setPCommandBuffers(&m_commandBuffers[imageIndex]);  // Command buffers(s) to execute in this batch (submission)
    submitInfo.setCommandBufferCount(1);                           // One command buffer
    submitInfo.setPNext(&deviceGroupSubmitInfo);

    // Submit to the graphics queue passing a wait fence
    m_queue.submit(submitInfo, m_waitFences[imageIndex]);

    // Presenting frame
    m_swapChain.present(m_queue);
  }


  //--------------------------------------------------------------------------------------------------
  // When the pipeline is set for using dynamic, this becomes useful
  //
  void setViewport(const vk::CommandBuffer& cmdBuf)
  {
    cmdBuf.setViewport(0, {vk::Viewport(0.0f, 0.0f, static_cast<float>(m_size.width), static_cast<float>(m_size.height), 0.0f, 1.0f)});
    cmdBuf.setScissor(0, {{{0, 0}, {m_size.width, m_size.height}}});
  }

  //--------------------------------------------------------------------------------------------------
  // Window callback when the it is resized
  // - Destroy allocated frames, then rebuild them with the new size
  // - Call onResize() of the derived class
  //
  virtual void onWindowResize(int w, int h)
  {
    if(w == 0 || h == 0)
    {
      return;
    }

    m_size.width  = w;
    m_size.height = h;

    // Update imgui and camera

    if(ImGui::GetCurrentContext() != nullptr)
    {
      auto& imgui_io       = ImGui::GetIO();
      imgui_io.DisplaySize = ImVec2(static_cast<float>(w), static_cast<float>(h));
    }
    CameraManip.setWindowSize(w, h);

    m_device.waitIdle();
    m_queue.waitIdle();

    m_swapChain.update(m_size.width, m_size.height, m_vsync);
    onResize(w, h);
    createDepthBuffer();
    createFrameBuffers();
  }


  //--------------------------------------------------------------------------------------------------
  // Window callback when the mouse move
  // - Handling ImGui and a default camera
  //
  virtual void onMouseMotion(int x, int y)
  {
    if(ImGui::GetCurrentContext() != nullptr && ImGui::GetIO().WantCaptureMouse)
    {
      return;
    }

    int ox, oy;  // Old mouse position
    CameraManip.getMousePosition(ox, oy);

    if(m_inputs.lmb || m_inputs.rmb || m_inputs.mmb)
    {
      CameraManip.mouseMove(x, y, m_inputs);
    }

    const int dx = x - ox, dy = y - oy;  // difference mouse position
    //---------------------------- LEFT
    if(m_inputs.lmb)
    {
      const float hval  = 2 * dx / static_cast<float>(m_size.width);
      const float vval  = 2 * dy / static_cast<float>(m_size.height);
      m_inertCamera.tau = s_cameraTau;
      m_inertCamera.rotateH(hval);
      m_inertCamera.rotateV(vval);
    }
    //---------------------------- MIDDLE
    if(m_inputs.mmb)
    {
      const float hval  = 2 * dx / static_cast<float>(m_size.width);
      const float vval  = 2 * dy / static_cast<float>(m_size.height);
      m_inertCamera.tau = s_cameraTau;
      m_inertCamera.rotateH(hval, true);
      m_inertCamera.rotateV(vval, true);
    }
    //---------------------------- RIGHT
    if(m_inputs.rmb)
    {
      const float hval  = 2 * dx / static_cast<float>(m_size.width);
      const float vval  = -2 * dy / static_cast<float>(m_size.height);
      m_inertCamera.tau = s_cameraTau;
      m_inertCamera.rotateH(hval, m_inputs.ctrl);
      m_inertCamera.move(vval, m_inputs.ctrl);
    }
  }

  //--------------------------------------------------------------------------------------------------
  // Window callback when a special key gets hit
  // - Handling ImGui and a default camera
  //
  virtual void onKeyboard(int key, int scancode, int action, int mods)
  {
    const bool capture = ImGui::GetCurrentContext() != nullptr && ImGui::GetIO().WantCaptureKeyboard;
    const bool pressed = action != GLFW_RELEASE;

    if(key == GLFW_KEY_LEFT_CONTROL)
    {
      m_inputs.ctrl = pressed;
    }
    if(key == GLFW_KEY_LEFT_SHIFT)
    {
      m_inputs.shift = pressed;
    }
    if(key == GLFW_KEY_LEFT_ALT)
    {
      m_inputs.alt = pressed;
    }

    // Skip key release, capture by ImGui except for F10 and ESC key
    if(action == GLFW_RELEASE || ((key == GLFW_KEY_F10 || key == GLFW_KEY_ESCAPE) ? false : capture))
    {
      return;
    }

    switch(key)
    {
      case GLFW_KEY_F10:
        m_show_gui = !m_show_gui;
        break;
      case GLFW_KEY_ESCAPE:
        glfwSetWindowShouldClose(m_window, 1);
        break;
      case GLFW_KEY_LEFT:
        m_inertCamera.tau = s_keyTau;
        m_inertCamera.rotateH(s_moveStep, m_inputs.ctrl);
        break;
      case GLFW_KEY_UP:
        m_inertCamera.tau = s_keyTau;
        m_inertCamera.rotateV(s_moveStep, m_inputs.ctrl);
        break;
      case GLFW_KEY_RIGHT:
        m_inertCamera.tau = s_keyTau;
        m_inertCamera.rotateH(-s_moveStep, m_inputs.ctrl);
        break;
      case GLFW_KEY_DOWN:
        m_inertCamera.tau = s_keyTau;
        m_inertCamera.rotateV(-s_moveStep, m_inputs.ctrl);
        break;
      case GLFW_KEY_PAGE_UP:
        m_inertCamera.tau = s_keyTau;
        m_inertCamera.move(s_moveStep, m_inputs.ctrl);
        break;
      case GLFW_KEY_PAGE_DOWN:
        m_inertCamera.tau = s_keyTau;
        m_inertCamera.move(-s_moveStep, m_inputs.ctrl);
        break;
      default:
        break;
    }
  }

  //--------------------------------------------------------------------------------------------------
  // Window callback when a key gets hit
  //
  virtual void onKeyboardChar(unsigned char key)
  {
    if(ImGui::GetCurrentContext() != nullptr && ImGui::GetIO().WantCaptureKeyboard)
    {
      return;
    }

    // Toggling vsync
    if(key == 'v')
    {
      m_vsync = !m_vsync;
      m_device.waitIdle();
      m_queue.waitIdle();
      m_swapChain.update(m_size.width, m_size.height, m_vsync);
      createFrameBuffers();
    }

    if(key == 'h' || key == '?')
      m_showHelp = !m_showHelp;
  }

  //--------------------------------------------------------------------------------------------------
  // Window callback when the mouse button is pressed
  // - Handling ImGui and a default camera
  //
  virtual void onMouseButton(int button, int action, int mods)
  {
    (void)mods;

    double x, y;
    glfwGetCursorPos(m_window, &x, &y);
    CameraManip.setMousePosition(static_cast<int>(x), static_cast<int>(y));

    m_inputs.lmb = (button == GLFW_MOUSE_BUTTON_LEFT) && (action == GLFW_PRESS);
    m_inputs.mmb = (button == GLFW_MOUSE_BUTTON_MIDDLE) && (action == GLFW_PRESS);
    m_inputs.rmb = (button == GLFW_MOUSE_BUTTON_RIGHT) && (action == GLFW_PRESS);
  }

  //--------------------------------------------------------------------------------------------------
  // Window callback when the mouse wheel is modified
  // - Handling ImGui and a default camera
  //
  virtual void onMouseWheel(int delta)
  {
    if(ImGui::GetCurrentContext() != nullptr && ImGui::GetIO().WantCaptureMouse)
    {
      return;
    }

    CameraManip.wheel(delta > 0 ? 1 : -1, m_inputs);

    m_inertCamera.tau = s_keyTau;
    m_inertCamera.move(delta > 0 ? s_moveStep : -s_moveStep, m_inputs.ctrl);
  }

  virtual void onFileDrop(const char* filename) {}

  //--------------------------------------------------------------------------------------------------
  // Initialization of the GUI
  // - Need to be call after the device creation
  //
  void initGUI(uint32_t subpassID = 0)
  {
    assert(m_renderPass && "Render Pass must be set");

    // UI
    ImGui::CreateContext();
    ImGuiIO& io    = ImGui::GetIO();
    io.IniFilename = nullptr;  // Avoiding the INI file
    io.LogFilename = nullptr;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    ImGuiH::setStyle();
    ImGuiH::setFonts();

    ImGui::InitVK(m_device, m_physicalDevice, m_queue, m_graphicsQueueIndex, m_renderPass, subpassID);
  }

  //--------------------------------------------------------------------------------------------------
  // Fit the camera to the Bounding box
  //
  void fitCamera(const nvmath::vec3f& boxMin, const nvmath::vec3f& boxMax, bool instantFit = true)
  {
    CameraManip.fit(boxMin, boxMax, instantFit, false, m_size.width / static_cast<float>(m_size.height));
  }

  // Return true if the window is minimized
  bool isMinimized(bool doSleeping = true)
  {
    int w, h;
    glfwGetWindowSize(m_window, &w, &h);
    bool minimized(w == 0 || h == 0);
    if(minimized && doSleeping)
    {
#ifdef _WIN32
      Sleep(50);
#else
      usleep(50);
#endif
    }
    return minimized;
  }

  void setTitle(const std::string& title) { glfwSetWindowTitle(m_window, title.c_str()); }

  // GLFW Callback setup
  void setupGlfwCallbacks(GLFWwindow* window)
  {
    m_window = window;
    glfwSetWindowUserPointer(window, this);
    glfwSetKeyCallback(window, &key_cb);
    glfwSetCharCallback(window, &char_cb);
    glfwSetCursorPosCallback(window, &cursorpos_cb);
    glfwSetMouseButtonCallback(window, &mousebutton_cb);
    glfwSetScrollCallback(window, &scroll_cb);
    glfwSetWindowSizeCallback(window, &windowsize_cb);
    glfwSetDropCallback(window, &drop_cb);
  }
  static void windowsize_cb(GLFWwindow* window, int w, int h)
  {
    auto app = reinterpret_cast<AppBase*>(glfwGetWindowUserPointer(window));
    app->onWindowResize(w, h);
  }
  static void mousebutton_cb(GLFWwindow* window, int button, int action, int mods)
  {
    auto app = reinterpret_cast<AppBase*>(glfwGetWindowUserPointer(window));
    app->onMouseButton(button, action, mods);
  }
  static void cursorpos_cb(GLFWwindow* window, double x, double y)
  {
    auto app = reinterpret_cast<AppBase*>(glfwGetWindowUserPointer(window));
    app->onMouseMotion(static_cast<int>(x), static_cast<int>(y));
  }
  static void scroll_cb(GLFWwindow* window, double x, double y)
  {
    auto app = reinterpret_cast<AppBase*>(glfwGetWindowUserPointer(window));
    app->onMouseWheel(static_cast<int>(y));
  }
  static void key_cb(GLFWwindow* window, int key, int scancode, int action, int mods)
  {
    auto app = reinterpret_cast<AppBase*>(glfwGetWindowUserPointer(window));
    app->onKeyboard(key, scancode, action, mods);
  }
  static void char_cb(GLFWwindow* window, unsigned int key)
  {
    auto app = reinterpret_cast<AppBase*>(glfwGetWindowUserPointer(window));
    app->onKeyboardChar(key);
  }
  static void drop_cb(GLFWwindow* window, int count, const char** paths)
  {
    auto app = reinterpret_cast<AppBase*>(glfwGetWindowUserPointer(window));
    int  i;
    for(i = 0; i < count; i++)
      app->onFileDrop(paths[i]);
  }
  // GLFW Callback end

  // Set if Nvlink will be used
  void useNvlink(bool useNvlink) { m_useNvlink = useNvlink; }

  //--------------------------------------------------------------------------------------------------
  // Getters
  vk::Instance                          getInstance() { return m_instance; }
  vk::Device                            getDevice() { return m_device; }
  vk::PhysicalDevice                    getPhysicalDevice() { return m_physicalDevice; }
  vk::Queue                             getQueue() { return m_queue; }
  uint32_t                              getQueueFamily() { return m_graphicsQueueIndex; }
  vk::CommandPool                       getCommandPool() { return m_cmdPool; }
  vk::RenderPass                        getRenderPass() { return m_renderPass; }
  vk::Extent2D                          getSize() { return m_size; }
  vk::PipelineCache                     getPipelineCache() { return m_pipelineCache; }
  vk::SurfaceKHR                        getSurface() { return m_surface; }
  const std::vector<vk::Framebuffer>&   getFramebuffers() { return m_framebuffers; }
  const std::vector<vk::CommandBuffer>& getCommandBuffers() { return m_commandBuffers; }
  uint32_t                              getCurFrame() const { return m_swapChain.getActiveImageIndex(); }
  vk::Format                            getColorFormat() const { return m_colorFormat; }
  vk::Format                            getDepthFormat() const { return m_depthFormat; }
  bool                                  showGui() { return m_show_gui; }

protected:
  uint32_t getMemoryType(uint32_t typeBits, const vk::MemoryPropertyFlags& properties) const
  {
    auto deviceMemoryProperties = m_physicalDevice.getMemoryProperties();
    for(uint32_t i = 0; i < deviceMemoryProperties.memoryTypeCount; i++)
    {
      if(((typeBits & (1 << i)) > 0) && (deviceMemoryProperties.memoryTypes[i].propertyFlags & properties) == properties)
      {
        return i;
      }
    }
    std::string err = "Unable to find memory type " + vk::to_string(properties);
    LOGE(err.c_str());
    assert(0);
    return ~0u;
  }

  // Showing help
  void uiDisplayHelp()
  {
    if(m_showHelp)
    {
      ImGui::BeginChild("Help", ImVec2(370, 120), true);
      ImGui::Text("%s", CameraManip.getHelp().c_str());
      ImGui::EndChild();
    }
  }

  //--------------------------------------------------------------------------------------------------

  // Vulkan low level
  vk::Instance       m_instance;
  vk::Device         m_device;
  vk::SurfaceKHR     m_surface;
  vk::PhysicalDevice m_physicalDevice;
  vk::Queue          m_queue;
  uint32_t           m_graphicsQueueIndex{VK_QUEUE_FAMILY_IGNORED};
  vk::CommandPool    m_cmdPool;

  // Drawing/Surface
  nvvk::SwapChain                m_swapChain;
  std::vector<vk::Framebuffer>   m_framebuffers;      // All framebuffers, correspond to the Swapchain
  std::vector<vk::CommandBuffer> m_commandBuffers;    // Command buffer per nb element in Swapchain
  std::vector<vk::Fence>         m_waitFences;        // Fences per nb element in Swapchain
  vk::Image                      m_depthImage;        // Depth/Stencil
  vk::DeviceMemory               m_depthMemory;       // Depth/Stencil
  vk::ImageView                  m_depthView;         // Depth/Stencil
  vk::RenderPass                 m_renderPass;        // Base render pass
  vk::Extent2D                   m_size{0, 0};        // Size of the window
  vk::PipelineCache              m_pipelineCache;     // Cache for pipeline/shaders
  bool                           m_vsync{false};      // Swapchain with vsync
  bool                           m_useNvlink{false};  // NVLINK usage
  GLFWwindow*                    m_window{nullptr};   // GLFW Window

  // Surface buffer formats
  vk::Format m_colorFormat{vk::Format::eB8G8R8A8Unorm};
  vk::Format m_depthFormat{vk::Format::eUndefined};

  // Two different camera manipulators
  nvh::CameraManipulator::Inputs m_inputs;       // Camera manipulator, like in Maya, 3dsmax, Softimage, ...
  InertiaCamera                  m_inertCamera;  // Camera Inertia

  // Other
  bool m_showHelp{false};  // Show help, pressing
  bool m_show_gui{true};
};  // namespace nvvk


}  // namespace nvvk
