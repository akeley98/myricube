/* Copyright (c) 2020, NVIDIA CORPORATION. All rights reserved.
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

#ifndef VKDO_FRAMEMANAGER_HPP_
#define VKDO_FRAMEMANAGER_HPP_

#include <array>
#include <cassert>
#include <deque>
#include <functional>
#include <stddef.h>
#include <stdexcept>
#include <stdint.h>

#include <vulkan/vulkan.h>
#include <nvvk/context_vk.hpp>
#include <nvvk/error_vk.hpp>
#include <nvvk/swapchain_vk.hpp>

namespace vkdo {

class FrameManager;

using FrameGarbageCallback = std::function<void(const FrameManager&)>;


// This class provides a higher-level abstraction around a Vulkan swap
// chain and command pool, and assists in setting up a simple
// double-buffer scheme for
//
// * synchronizing resources between the host and the presentation queue,
//
// * managing freeing resources used for drawing a frame, and
//
// * allocating and auto-freeing single-use command buffers submitted
//   to the presentation queue.
//
// This class is not designed to handle resources used or shared with
// queues besides the presentation queue used to initialize this
// class; however, this class is cooperative in that it never blocks
// any queue other than its own presentation queue (i.e. no
// vkDeviceWaitIdle).
//
// To use this class, instantiate it with the nvvk::Context, surface
// to draw to, and the surface's width and height. This adapts and
// uses the GCT queue of the context for presentation (an alternate
// constructor is provided that doesn't use nvvk::Context and allows a
// manual queue choice).
//
// Then, make 2 copies of every resource you want to synchronize; call
// one set the even set and one the odd set. This class ensures at
// most 2 consecutive frames are in-flight at once, so, all
// odd-numbered frames safely re-use the same odd resource set (and
// same for the even set).
//
// A frame starts with a call to beginFrame, which returns a primary
// command buffer (among other things) and ends with endFrame, which
// submits said command buffer to the presentation queue. Within a
// beginFrame/endFrame pair, you may:
//
// * Use evenOdd() to select between resources in the even set or odd
//   set, as appropriate for this frame.
//
//   Within the beginFrame/endFrame pair, only 1 frame is in flight
//   (of opposite parity as the currently recorded frame), so the
//   returned object can safely be manipulated by the host.
//
// * Use addFrameGarbage() to schedule stuff for destruction, which
//   will occur only once the currently-recorded frame finishes
//   execution (i.e. 2 frames from now, or when the class destructor runs).
//
//   TODO consider writing easier/more efficient versions for common resources.
//
// Finally, as an exception to double-buffering, the class destructor,
// and beginFrame() IF it recreates the swap chain, block the
// presentation queue entirely. This simplifies tasks like cleanup and
// recreating frame buffers.
class FrameManager
{
    // Number of frames started since construction (so first frame is frame 1).
    uint64_t m_frameNumber = 0;

    // The Vulkan instance + device this FrameManager is constructed for.
    VkInstance m_instance;
    VkPhysicalDevice m_physicalDevice;
    VkDevice m_device;

    // Must be set to a queue that supports graphics, compute,
    // transfer, and presenting to the swap chain (constructor from
    // nvvk::Context takes care of this). Will be used for presenting
    // to the swap chain (i.e. the window).
    VkQueue m_presentQueue;

    // Objects above are borrowed; we own the ones below and create +
    // destroy them.

    // Command Pool for allocating one-time command buffers suitable
    // for the above queue.
    VkCommandPool m_oneTimeCommandPool;

    // The command buffer given to the user by
    // beginFrame. VK_NULL_HANDLE when not between a
    // beginFrame/endFrame pair.
    VkCommandBuffer m_userCommandBuffer = VK_NULL_HANDLE;

    // Abstracts away most swap chain stuff. We store the width and
    // height of the actual swap chain image (which may differ from
    // requested).
    nvvk::SwapChain m_swapChain;
    uint32_t m_width, m_height;

    // m_frameFences[0] is signalled when an even-numbered frame is
    // finished (by the device), m_frameFences[1] for odd.
    //
    // These fences must always either be signalled, or scheduled to
    // be signalled, except that in between beginFrame and endFrame,
    // the fence corresponding to that frame will not be signalled.
    // This means it's always safe to unconditionally wait on a fence
    // before starting a new frame.
    std::array<VkFence, 2> m_frameFences;

    // Lists of garbage to be destroyed.  Stuff is pushed onto
    // m_garbageLists[0] or [1] depending on the parity of the frame.
    // The list will be destroyed after we wait on the fence for that
    // frame, in front-to-back order.
    std::array<std::deque<FrameGarbageCallback>, 2> m_garbageLists;

    // In case anyone wants to write their own malloc.
    VkAllocationCallbacks* pAllocator = nullptr;

    void destroyGarbageList(std::deque<FrameGarbageCallback>& garbage) noexcept
    {
        for (FrameGarbageCallback& f : garbage) f(*this);
        garbage.clear();
    }

    static constexpr uint64_t forever = ~uint64_t(0);

  public:
    // Manual constructor: you pass in
    //
    // The vulkan instance, device,
    // and physical device you want to use
    //
    // A surface to render to and its dimensions
    //
    // A queue (plus its queue family index) that MUST be capable of
    // drawing to the surface, and graphics, compute, and transfer operations.
    FrameManager(
        VkInstance instance,
        VkPhysicalDevice physicalDevice, VkDevice device,
        VkSurfaceKHR surface, uint32_t width, uint32_t height,
        VkQueue queue, uint32_t queueFamilyIndex) :
    m_instance(instance),
    m_physicalDevice(physicalDevice),
    m_device(device),
    m_presentQueue(queue),
    m_width(width), m_height(height)
    {
        this->constructor(surface, queueFamilyIndex);
    }

    // Constructor from nvvk::Context. Sets the GCT queue of the
    // context to one usable by the surface, then steals it for
    // ourselves. You stil have to provide the surface to render to
    // and its dimensions.
    //
    // (Consider glfwCreateWindowSurface and glfwGetFramebufferSize).
    FrameManager(
        nvvk::Context& ctx,
        VkSurfaceKHR surface, uint32_t width, uint32_t height) :
    m_instance(ctx.m_instance),
    m_physicalDevice(ctx.m_physicalDevice),
    m_device(ctx.m_device),
    m_presentQueue((ctx.setGCTQueueWithPresent(surface), ctx.m_queueGCT.queue)),
    m_width(width), m_height(height)
    {
        // ^^^ Comma trick above keeps m_presentQueue as const.
        uint32_t queueFamilyIndex = ctx.m_queueGCT.familyIndex;
        this->constructor(surface, queueFamilyIndex);
    }

  private:
    void constructor(VkSurfaceKHR surface, uint32_t queueFamilyIndex)
    {
        // Set up the command pool.
        VkCommandPoolCreateInfo commandPoolArgs{
            VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO,
            nullptr, 0, queueFamilyIndex };
        NVVK_CHECK(vkCreateCommandPool(
            m_device, &commandPoolArgs, pAllocator, &m_oneTimeCommandPool));

        // Swap chain needs to be manually initialized.
        m_swapChain.init(
            m_device, m_physicalDevice,
            m_presentQueue, queueFamilyIndex,
            surface,
            VK_FORMAT_B8G8R8A8_UNORM); // Document format?

        // This is needed as we promised we wouldn't block the whole
        // device (thread safety). No vkDeviceWaitIdle.
        m_swapChain.setWaitQueue(m_presentQueue);

        m_swapChain.update(m_width, m_height);

        // Initialize fences in signalled state as specified.
        VkFenceCreateInfo fenceArgs {
            VK_STRUCTURE_TYPE_FENCE_CREATE_INFO,
            nullptr,
            VK_FENCE_CREATE_SIGNALED_BIT };
        for (int i = 0; i < 2; ++i) {
            VkFence* pFence = &m_frameFences[i];
            NVVK_CHECK(vkCreateFence(m_device, &fenceArgs, pAllocator, pFence));
        }
    }

  public:
    FrameManager(FrameManager&&) = delete;

    ~FrameManager()
    {
        assert(!inBeginEndPair());

        // Wait and destroy fences.
        vkWaitForFences(m_device, 2, m_frameFences.data(), VK_TRUE, forever);
        vkDestroyFence(m_device, m_frameFences[0], pAllocator);
        vkDestroyFence(m_device, m_frameFences[1], pAllocator);

        // Now that the fences are gone, we're safe to destroy
        // everything.  Always destroy the older frame's stuff first,
        // to reduce unpredictability. Keep in mind beginFrame, not
        // endFrame, bumps m_frameNumber.
        destroyGarbageList(m_garbageLists[1 ^ (m_frameNumber & 1)]);
        destroyGarbageList(m_garbageLists[m_frameNumber & 1]);

        // Finally destroy the other stuff we own.
        m_swapChain.deinit();
        vkDestroyCommandPool(m_device, m_oneTimeCommandPool, pAllocator);
    }

    // Select one reference arg or the other depending on the parity
    // of the current frame. The one NOT returned must not be
    // modified or destroyed by the host, as it can be in-use by the
    // device.
    template <typename T>
    T&& evenOdd(T&& useOnEven, T&& useOnOdd) const
    {
        assert(inBeginEndPair());
        return std::forward<T>(m_frameNumber & 1 ? useOnEven : useOnOdd);
    }

    template <typename T, size_t Two>
    T& evenOdd(T array[Two]) const
    {
        static_assert(Two == 2, "Must be size 2 array.");
        return array[evenOdd()];
    }

    template <typename Container>
    auto& evenOdd(Container& container) const
    {
        assert(container.size() == 2);
        return container[evenOdd()];
    }

    // Return a raw 0 or 1 value depending on whether the frame number
    // is even or odd.
    int evenOdd() const
    {
        assert(inBeginEndPair());
        return int(m_frameNumber & 1);
    }

    // Return whether we're in-between a beginFrame/endFrame
    bool inBeginEndPair() const
    {
        return m_userCommandBuffer != VK_NULL_HANDLE;
    }

    // Allocate a primary command buffer and start its recording. This
    // command buffer is suitable for submitting to m_presentQueue
    // exactly once.
    VkCommandBuffer recordOneTimeCommandBuffer()
    {
        VkCommandBufferAllocateInfo allocInfo {
            VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,
            nullptr,
            m_oneTimeCommandPool,
            VK_COMMAND_BUFFER_LEVEL_PRIMARY,
            1 };
        VkCommandBuffer cmdBuffer;
        NVVK_CHECK(vkAllocateCommandBuffers(m_device, &allocInfo, &cmdBuffer));

        VkCommandBufferBeginInfo beginInfo {
            VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,
            nullptr,
            VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT,
            nullptr };
        NVVK_CHECK(vkBeginCommandBuffer(cmdBuffer, &beginInfo));

        return cmdBuffer;
    }

    // TODO secondary command buffers.


    // Start a new frame.
    //
    // * Increment the current frame number.
    //
    // * Wait for the frame 2 frames ago to finish, and clean up
    //   garbage registered that frame. Now the even (odd) resource
    //   set is safe to modify if the current frame number is even (odd).
    //
    // * Allocate a one-time command buffer for the presentation queue
    //   and return it through *pCmdBuffer.
    //
    // * Acquire a new swap chain image, and return it through
    //   *pSwapChainImage.
    //
    // * Read the requested swap chain width/height from *pWidth and
    //   *pHeight (this should be window framebuffer size; I use
    //   glfwGetFramebufferSize), and overwrite *pWidth and *pHeight
    //   with the actual swap image size (which may be different!!!)
    //
    // * Re-create the swap chain if needed (typically due to screen
    //   size change, but could happen unexpectedly); if so,
    //   vkQueueWaitIdle is called for the presentation queue. Return
    //   through *pSwapChainRecreated whether this re-create happened.
    //
    // Except for pSwapChainRecreated, all pointers must be non-nullptr.
    void beginFrame(
        VkCommandBuffer* pCmdBuffer,
        nvvk::SwapChainImage* pSwapChainImage,
        uint32_t* pWidth, uint32_t* pHeight,
        bool* pSwapChainRecreated=nullptr) noexcept
    {
        // Increment frame counter.
        assert(!inBeginEndPair());
        ++m_frameNumber; // Exception safety? set noexcept for now

        // Record a new command buffer for this frame.
        assert(pCmdBuffer != nullptr);
        *pCmdBuffer = m_userCommandBuffer = recordOneTimeCommandBuffer();

        // Wait for the frame 2 frames ago to finish, then clean up
        // its garbage. Need to wait before asking for swap image.
        VkFence frameFence = evenOdd(m_frameFences);
        NVVK_CHECK(vkWaitForFences(m_device, 1, &frameFence, VK_TRUE, forever));
        NVVK_CHECK(vkResetFences(m_device, 1, &frameFence));
        destroyGarbageList(evenOdd(m_garbageLists));

        // Get the next swap chain image.
        assert(pSwapChainImage != nullptr);
        assert(pWidth != nullptr);
        assert(pHeight != nullptr);
        m_swapChain.acquireAutoResize(
            *pWidth, *pHeight, pSwapChainRecreated, pSwapChainImage);

        // Return the actual swap chain image size.
        VkExtent2D extent = m_swapChain.getExtent();
        *pWidth  = extent.width;
        m_width  = extent.width;
        *pHeight = extent.height;
        m_height = extent.height;
    }

    // Record a command for transitioning the layout of the current
    // swap chain image from the given oldLayout to
    // VK_IMAGE_LAYOUT_PRESENT_SRC_KHR. This also defines a memory
    // barrier operation ensuring all writes (done on this queue) to
    // the swap chain image finish before layout transition (unless
    // you pass accessFlags manually).
    //
    // Requires that swap chain image is owned by m_presentQueue (if
    // you don't use multiple queues, it is).
    void cmdSwapChainImageFixLayout(
        VkCommandBuffer cmdBuf,
        VkImageLayout oldLayout,
        VkAccessFlagBits accessFlags = VK_ACCESS_MEMORY_WRITE_BIT,
        VkPipelineStageFlags stageFlags = VK_PIPELINE_STAGE_ALL_COMMANDS_BIT)
    {
        VkImageMemoryBarrier imageLayoutBarrier {
            VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER,
            nullptr,
            accessFlags,
            VK_ACCESS_MEMORY_READ_BIT,
            oldLayout,
            VK_IMAGE_LAYOUT_PRESENT_SRC_KHR,
            VK_QUEUE_FAMILY_IGNORED,
            VK_QUEUE_FAMILY_IGNORED,
            m_swapChain.getActiveImage(),
            { VK_IMAGE_ASPECT_COLOR_BIT,
              0, VK_REMAINING_MIP_LEVELS,
              0, VK_REMAINING_ARRAY_LAYERS } };

        vkCmdPipelineBarrier(
            cmdBuf,
            stageFlags,
            VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT,
            0, 0, nullptr,
            0, nullptr,
            1, &imageLayoutBarrier);
    }

    // End the started frame. Submits the command buffer returned by
    // beginFrame to the presentation queue, schedules it to be freed
    // later, and present the swap chain image acquired in beginFrame.
    //
    // userCommandBuffer must be the same command buffer that the last
    // call to beginFrame returned. You must NOT free this command
    // buffer manually.
    void endFrame(VkCommandBuffer userCommandBuffer)
    {
        // Finish recording the command buffer, which must be the
        // one beginFrame gave out.
        assert(userCommandBuffer == m_userCommandBuffer);
        assert(inBeginEndPair());
        NVVK_CHECK(vkEndCommandBuffer(userCommandBuffer));

        // Submit command buffer to queue, and signal the correct
        // fence for this frame's parity.
        VkFence frameFence = evenOdd(m_frameFences);
        VkSemaphore waitSem = m_swapChain.getActiveReadSemaphore();
        VkSemaphore signalSem = m_swapChain.getActiveWrittenSemaphore();

        VkPipelineStageFlags semStageMask = VK_PIPELINE_STAGE_ALL_COMMANDS_BIT;
        VkSubmitInfo submitInfo {
            VK_STRUCTURE_TYPE_SUBMIT_INFO,
            nullptr,
            1, &waitSem, &semStageMask,
            1, &userCommandBuffer,
            1, &signalSem };
        NVVK_CHECK(vkQueueSubmit(m_presentQueue, 1, &submitInfo, frameFence));

        // Schedule the CommandBuffer for later destruction.
        addFrameGarbage([userCommandBuffer] (const FrameManager& self) {
            vkFreeCommandBuffers(
                self.getDevice(),
                self.getCommandPool(),
                1, &userCommandBuffer);
        });
        m_userCommandBuffer = VK_NULL_HANDLE;

        // Present the drawn image.
        m_swapChain.present();
    }

    // Schedule this callback to be called when the current frame is
    // finished on the device. This can be used to do arbitrary work,
    // but the design case was for dealing with garbage (single use
    // command buffers, etc.)
    // Must be called only in beginFrame/endFrame pairs.
    //
    // Garbage callbacks are called in reverse-order of their
    // registration, matching typically expected behavior (C++
    // destructors, atexit, etc.)
    //
    // If you want to write a fancy C++ callback, use the next
    // function, which just takes an object and a function pointer.
    void addFrameGarbage(FrameGarbageCallback garbage)
    {
        evenOdd(m_garbageLists).push_front(std::move(garbage));
    }

    template <typename T> void addFrameGarbage(
        T obj, void (*destroyer) (T, const FrameManager&))
    {
        addFrameGarbage([obj, destroyer] (const FrameManager& self) {
            destroyer(obj, self); } );
    }

    // Like the above two functions, but callbacks are called in the
    // order they're registered.
    void addFrameGarbageLast(FrameGarbageCallback garbage)
    {
        evenOdd(m_garbageLists).push_back(std::move(garbage));
    }

    template <typename T> void addFrameGarbageLast(
        T obj, void (*destroyer) (T, const FrameManager&))
    {
        addFrameGarbageLast([obj, destroyer] (const FrameManager& self) {
            destroyer(obj, self); } );
    }

    // Dumb getters.
    uint64_t getFrameNumber() const { return m_frameNumber; }

    VkInstance getInstance() const { return m_instance; }
    VkPhysicalDevice getPhysicalDevice() const { return m_physicalDevice; }
    VkDevice getDevice() const { return m_device; }

    VkQueue getQueue() const { return m_presentQueue; }
    VkCommandPool getCommandPool() const { return m_oneTimeCommandPool; }

    void getWidthHeight(uint32_t* pWidth, uint32_t* pHeight) const
    {
        if (pWidth) *pWidth = m_width;
        if (pHeight) *pHeight = m_height;
    }

    nvvk::SwapChain& getSwapChain(
        uint32_t* pWidth=nullptr, uint32_t* pHeight=nullptr)
    {
        getWidthHeight(pWidth, pHeight);
        return m_swapChain;
    }
};

} // end namespace

#endif /* end include guard */
