// Vulkan implementation of RendererLogic.

#ifdef WIN32
    #define VK_USE_PLATFORM_WIN32_KHR
#else
    #define VK_USE_PLATFORM_XCB_KHR
    #define VK_USE_PLATFORM_XLIB_KHR
#endif

#define GLFW_INCLUDE_VULKAN
#include "GLFW/glfw3.h"

#include "RendererLogic.hh"

#include <cassert>
#include <vulkan/vulkan.h>

#include "nvvk/context_vk.hpp"
#include "FrameManager.hpp"

using namespace myricube;
using namespace akeley;

namespace {

// Handle for GPU resources for the mesh of one chunk group.
struct MeshEntry
{

};

struct MeshStaging
{

};

// aabb_array[z][y][x] is the minimal AABB containing the visible
// voxels of ChunkGroup::chunk_array[z][y][x].
//
// If you change this, be aware that sizeof is used on this to
// know the amount of bytes to allocate for the GPU's VBO.
struct AABBs
{
    PackedAABB aabb_array[edge_chunks][edge_chunks][edge_chunks];
};

// All voxels within a chunk group share a single 3D voxel array on
// the GPU (stored as a 3D texture), as well as one VBO used to store
// the array of AABB for the chunks within the group. This is the
// handle for the GPU data needed to raycast one chunk group.
struct RaycastEntry
{

};

// Staging buffer for the async cache of raycast chunk groups.
struct RaycastStaging
{

};

// nvvk::Context with proper constructor and destructor.
struct ScopedContext : nvvk::Context
{
    ScopedContext()
    {
        // Add the swapchain extensions.
        nvvk::ContextCreateInfo deviceInfo;
        deviceInfo.apiMajor = 1;
        deviceInfo.apiMinor = 1;
        deviceInfo.addInstanceExtension(VK_KHR_SURFACE_EXTENSION_NAME);
#ifdef WIN32
        deviceInfo.addInstanceExtension(VK_KHR_WIN32_SURFACE_EXTENSION_NAME);
#else
        deviceInfo.addInstanceExtension(VK_KHR_XLIB_SURFACE_EXTENSION_NAME);
        deviceInfo.addInstanceExtension(VK_KHR_XCB_SURFACE_EXTENSION_NAME);
#endif
        deviceInfo.addDeviceExtension(VK_KHR_SWAPCHAIN_EXTENSION_NAME);

        // Initialize the nvvk context.
        init(deviceInfo);
    }

    ScopedContext(ScopedContext&&) = delete;

    ~ScopedContext()
    {
        deinit();
    }
};

// Thingie for managing GLFW's VkSurfaceKHR.
struct Surface
{
    // Managed by us.
    VkSurfaceKHR surface;

    // Borrowed pointers.
    GLFWwindow* p_window;
    VkInstance instance;

    int initial_width;
    int initial_height;

    Surface(VkInstance inst_, GLFWwindow* window)
    {
        p_window = window;
        instance = inst_;
        NVVK_CHECK(glfwCreateWindowSurface(inst_, window, nullptr, &surface));
        glfwGetFramebufferSize(window, &initial_width, &initial_height);
    }

    Surface(Surface&&) = delete;

    ~Surface()
    {
        vkDestroySurfaceKHR(instance, surface, nullptr);
        surface = VK_NULL_HANDLE;
        instance = VK_NULL_HANDLE;
    }
};

} // end anonymous namespace.



namespace myricube {

struct RendererVk :
    RendererLogic<MeshEntry, MeshStaging, RaycastEntry, RaycastStaging>
{
    ScopedContext ctx;
    Surface glfw_surface;
    FrameManager frame_manager;

    // set in begin_frame.
    int width, height;

    VkCommandBuffer frame_cmd_buffer = VK_NULL_HANDLE;
    nvvk::SwapChainImage current_swap_image;

    RendererVk(RenderThread* thread, RenderArgs args) :
        RendererLogic<MeshEntry, MeshStaging, RaycastEntry, RaycastStaging>(
            thread,
            args),
        glfw_surface(
            ctx.m_instance,
            args.p_window->get_glfw_window()),
        frame_manager(
            ctx,
            glfw_surface.surface,
            glfw_surface.initial_width,
            glfw_surface.initial_height)
    {

    }

    void begin_frame() override
    {
        glfwGetFramebufferSize(glfw_surface.p_window, &width, &height);
        frame_manager.beginFrame(
            &frame_cmd_buffer,
            &current_swap_image,
            width, height);

        VkImageSubresourceRange imageRange {
            VK_IMAGE_ASPECT_COLOR_BIT,
            0, VK_REMAINING_MIP_LEVELS,
            0, VK_REMAINING_ARRAY_LAYERS };

        // Command swap chain image transition to dst optimal layout.
        VkImageMemoryBarrier transferLayoutBarrier {
            VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER,
            nullptr,
            0,
            VK_ACCESS_MEMORY_WRITE_BIT,
            VK_IMAGE_LAYOUT_UNDEFINED,
            VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
            VK_QUEUE_FAMILY_IGNORED,
            VK_QUEUE_FAMILY_IGNORED,
            current_swap_image.image,
            imageRange };

        vkCmdPipelineBarrier(
            frame_cmd_buffer,
            VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            0, 0, nullptr,
            0, nullptr,
            1, &transferLayoutBarrier);

        // Clear the swap chain image to one color.
        uint64_t ns = uint64_t(glfwGetTime() * 1e9);

        VkClearColorValue color{};
        double nanoTau = 6.283185307179587e-9;
        float red = 0.5 + 0.5 * cos((ns % 1'000'000'000) * nanoTau);
        float green = 0.5 + 0.5 * cos(2.0 + (ns % 1'000'000'000) * nanoTau);
        float blue = 1.0f - red - green;
        color.float32[0] = red;
        color.float32[1] = green;
        color.float32[2] = blue;
        color.float32[3] = 1.0f;
        vkCmdClearColorImage(
            frame_cmd_buffer,
            current_swap_image.image,
            VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
            &color,
            1, &imageRange);

        // Transition image to layout suitable for display.
        frame_manager.cmdSwapChainImageFixLayout(
            frame_cmd_buffer,
            VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
            VK_ACCESS_TRANSFER_WRITE_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT);
    }


    /* MESH RENDERING IMPLEMENTATION */

    // Draw the list of chunk groups with the instanced-voxel
    // rendering method.
    void draw_mesh_entries(
        const std::vector<std::pair<MeshEntry*, glm::ivec3>>& entries) override
    {

    }

    // Convert chunk group's chunks to meshes and pack into the
    // memory-mapped array of chunks.
    void worker_stage(
        MeshStaging* staging, const BinChunkGroup* group_ptr) override
    {

    }

    // Since MeshEntry and MeshStaging are actually the same,
    // literally swap-in the staged buffer into the cache (with the
    // evicted cache entry, which is no longer used for rendering,
    // used as the new staging buffer).
    bool swap_in(
        MeshStaging* staging,
        std::unique_ptr<MeshEntry>* p_uptr_entry) override
    {
        return true;
    }



    /* RAYCAST RENDERER IMPLEMENTATION */

    // Big-picture of data: I'm storing voxel data for chunk groups as
    // 3D VkImages. These textures live in RaycastEntry inside the
    // main cache of AsyncCache, hidden in RendererLogic.
    //
    // To fill these textures, worker threads in AsyncCache fill the
    // staging buffers (memory-mapped SSBO) with voxel data loaded
    // from disk. This is separate from the main cache.
    //
    // Once the worker thread marks a staging buffer as fully loaded,
    // it can be swapped into the main cache. This is done in two steps:
    // the worker thread enqueues commands to the transfer-only queue
    // to copy from the staging buffer to VkImage, then the main thread
    // checks a VkFence to see if this is done and swaps the VkImage
    // into the main cache when ready.
    //
    // Both threads also need to schedule a pipeline barrier in order
    // to do the image layout transition and queue ownership transfer.
    // TODO implement

    // Draw the given list of chunk groups with the raycast-AABB
    // rendering method.
    void draw_raycast_entries(
        const std::vector<std::pair<RaycastEntry*, glm::ivec3>>& entries)
    override
    {

    }

    void worker_stage(
        RaycastStaging* stage, const BinChunkGroup* group_ptr) override
    {

    }

    bool swap_in(RaycastStaging* staging,
                 std::unique_ptr<RaycastEntry>* p_uptr_entry) override
    {
        return true;
    }



    void end_frame() override
    {
        frame_manager.endFrame(frame_cmd_buffer);
    }

    void wait_idle() override
    {
        vkDeviceWaitIdle(ctx.m_device);
    }
};

std::shared_ptr<RendererBase> RendererVk_Factory(
    RenderThread* thread,
    RenderArgs args)
{
    return std::shared_ptr<RendererBase>(new RendererVk(thread, args));
}


} // end namespace myricube
