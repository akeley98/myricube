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
#include <fstream>
#include <vulkan/vulkan.h>

#include "nvvk/context_vk.hpp"
#include "FrameManager.hpp"

using namespace myricube;
using namespace akeley;

namespace {

// These should be widely supported, but I can detect support if
// needed later.
constexpr auto swap_chain_image_format = VK_FORMAT_B8G8R8A8_SRGB;
constexpr auto depth_format = VK_FORMAT_D24_UNORM_S8_UINT;

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
struct ScopedSurface
{
    // Managed by us.
    VkSurfaceKHR surface;

    // Borrowed pointers.
    GLFWwindow* p_window;
    VkInstance instance;

    int initial_width;
    int initial_height;

    ScopedSurface(VkInstance inst_, GLFWwindow* window)
    {
        p_window = window;
        instance = inst_;
        NVVK_CHECK(glfwCreateWindowSurface(inst_, window, nullptr, &surface));
        glfwGetFramebufferSize(window, &initial_width, &initial_height);
    }

    ScopedSurface(ScopedSurface&&) = delete;

    ~ScopedSurface()
    {
        vkDestroySurfaceKHR(instance, surface, nullptr);
        surface = VK_NULL_HANDLE;
        instance = VK_NULL_HANDLE;
    }

    operator VkSurfaceKHR() const { return surface; }
};



// Thingie for managing the simple one-subpass, depth buffer
// VkRenderPass. Mostly copied from vulkan-tutorial.com.
struct ScopedRenderPass
{
    // Managed by us.
    VkRenderPass render_pass;

    // Borrowed pointers.
    VkDevice device;

    ScopedRenderPass(VkDevice dev_)
    {
        device = dev_;

        VkAttachmentDescription colorAttachment{};
        colorAttachment.format = swap_chain_image_format;
        colorAttachment.samples = VK_SAMPLE_COUNT_1_BIT;
        colorAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
        colorAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
        colorAttachment.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
        colorAttachment.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
        colorAttachment.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        colorAttachment.finalLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;

        VkAttachmentDescription depthAttachment{};
        depthAttachment.format = depth_format;
        depthAttachment.samples = VK_SAMPLE_COUNT_1_BIT;
        depthAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
        depthAttachment.storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
        depthAttachment.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
        depthAttachment.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
        depthAttachment.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        depthAttachment.finalLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;

        VkAttachmentReference colorAttachmentRef{};
        colorAttachmentRef.attachment = 0;
        colorAttachmentRef.layout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;

        VkAttachmentReference depthAttachmentRef{};
        depthAttachmentRef.attachment = 1;
        depthAttachmentRef.layout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;

        VkSubpassDescription subpass{};
        subpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
        subpass.colorAttachmentCount = 1;
        subpass.pColorAttachments = &colorAttachmentRef;
        subpass.pDepthStencilAttachment = &depthAttachmentRef;

        VkSubpassDependency dependency{};
        dependency.srcSubpass = VK_SUBPASS_EXTERNAL;
        dependency.dstSubpass = 0;
        dependency.srcStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        dependency.srcAccessMask = 0;
        dependency.dstStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        dependency.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;

        std::array<VkAttachmentDescription, 2> attachments = {colorAttachment, depthAttachment};
        VkRenderPassCreateInfo renderPassInfo{};
        renderPassInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO;
        renderPassInfo.attachmentCount = static_cast<uint32_t>(attachments.size());
        renderPassInfo.pAttachments = attachments.data();
        renderPassInfo.subpassCount = 1;
        renderPassInfo.pSubpasses = &subpass;
        renderPassInfo.dependencyCount = 1;
        renderPassInfo.pDependencies = &dependency;

        NVVK_CHECK(
            vkCreateRenderPass(device, &renderPassInfo, nullptr, &render_pass));
    }

    ScopedRenderPass(ScopedRenderPass&&) = delete;

    ~ScopedRenderPass()
    {
        vkDestroyRenderPass(device, render_pass, nullptr);
    }

    operator VkRenderPass() const
    {
        return render_pass;
    }
};



// Copied from vulkan-tutorial.com
std::vector<char> readFile(const std::string& filename)
{
    std::ifstream file(filename, std::ios::ate | std::ios::binary);

    if (!file.is_open()) {
        throw std::runtime_error("failed to open file: " + filename);
    }

    size_t fileSize = (size_t) file.tellg();
    std::vector<char> buffer(fileSize);

    file.seekg(0);
    file.read(buffer.data(), fileSize);

    file.close();

    return buffer;
}

VkShaderModule createShaderModule(VkDevice device, const std::vector<char>& code)
{
    VkShaderModuleCreateInfo createInfo{};
    createInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
    createInfo.codeSize = code.size();
    createInfo.pCode = reinterpret_cast<const uint32_t*>(code.data());

    VkShaderModule shaderModule;
    if (vkCreateShaderModule(device, &createInfo, nullptr, &shaderModule) != VK_SUCCESS) {
        throw std::runtime_error("failed to create shader module!");
    }

    return shaderModule;
}



struct PushConstant
{
    glm::mat4 mvp;
    glm::vec4 eye_relative_group_origin; // w unused.
    int32_t flags;
    int32_t far_plane_squared;
    int32_t raycast_thresh_squared;
};
static_assert(sizeof(PushConstant) <= 128);



// Pipeline for mesh rendering. All we really have to do is load the
// vert/frag shader and hook up the instanced voxels array.
struct MeshPipeline
{
    // We manage these.
    VkPipeline pipeline;
    VkPipelineLayout layout;

    // Borrowed pointer.
    VkDevice device;

    const uint32_t width, height;

    MeshPipeline(
        VkDevice dev_,
        VkRenderPass render_pass,
        uint32_t width_, uint32_t height_) :
    device(dev_),
    width(width_),
    height(height_)
    {
        // Boilerplate from vulkan-tutorial.com with some modifications.
        VkPushConstantRange pushConstantRange{};
        pushConstantRange.stageFlags = VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT;
        pushConstantRange.offset = 0;
        pushConstantRange.size = sizeof(PushConstant);

        VkPipelineLayoutCreateInfo pipelineLayoutInfo{};
        pipelineLayoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        pipelineLayoutInfo.setLayoutCount = 0;
        pipelineLayoutInfo.pSetLayouts = nullptr;
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushConstantRange;

        NVVK_CHECK(vkCreatePipelineLayout(device, &pipelineLayoutInfo, nullptr, &layout));

        auto vertShaderCode = readFile(expand_filename("vk/mesh.vert.spv"));
        auto fragShaderCode = readFile(expand_filename("vk/mesh.frag.spv"));

        VkShaderModule vertShaderModule = createShaderModule(device, vertShaderCode);
        VkShaderModule fragShaderModule = createShaderModule(device, fragShaderCode);

        VkPipelineShaderStageCreateInfo vertShaderStageInfo{};
        vertShaderStageInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        vertShaderStageInfo.stage = VK_SHADER_STAGE_VERTEX_BIT;
        vertShaderStageInfo.module = vertShaderModule;
        vertShaderStageInfo.pName = "main";

        VkPipelineShaderStageCreateInfo fragShaderStageInfo{};
        fragShaderStageInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        fragShaderStageInfo.stage = VK_SHADER_STAGE_FRAGMENT_BIT;
        fragShaderStageInfo.module = fragShaderModule;
        fragShaderStageInfo.pName = "main";

        VkPipelineShaderStageCreateInfo shaderStages[] = {vertShaderStageInfo, fragShaderStageInfo};

        VkPipelineVertexInputStateCreateInfo vertexInputInfo{};
        vertexInputInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO;

        VkVertexInputBindingDescription bindingDescription {
            0, // binding
            sizeof(MeshVoxelVertex), // stride
            VK_VERTEX_INPUT_RATE_INSTANCE, // Voxels array is instanced.
        };
        std::array<VkVertexInputAttributeDescription, 2> attributeDescriptions {
            VkVertexInputAttributeDescription {
              0, // location
              0, // binding
              VK_FORMAT_R32_UINT,
              offsetof(MeshVoxelVertex, packed_residue_face_bits) },
            VkVertexInputAttributeDescription {
              1, // location
              0, // binding
              VK_FORMAT_R32_UINT,
              offsetof(MeshVoxelVertex, packed_color) },
        };

        vertexInputInfo.vertexBindingDescriptionCount = 1;
        vertexInputInfo.vertexAttributeDescriptionCount = static_cast<uint32_t>(attributeDescriptions.size());
        vertexInputInfo.pVertexBindingDescriptions = &bindingDescription;
        vertexInputInfo.pVertexAttributeDescriptions = attributeDescriptions.data();

        VkPipelineInputAssemblyStateCreateInfo inputAssembly{};
        inputAssembly.sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO;
        inputAssembly.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
        inputAssembly.primitiveRestartEnable = VK_FALSE;

        VkViewport viewport{};
        viewport.x = 0.0f;
        viewport.y = 0.0f;
        viewport.width = width;
        viewport.height = height;
        viewport.minDepth = 0.0f;
        viewport.maxDepth = 1.0f;

        VkRect2D scissor{};
        scissor.offset = {0, 0};
        scissor.extent = VkExtent2D{ width, height };

        VkPipelineViewportStateCreateInfo viewportState{};
        viewportState.sType = VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO;
        viewportState.viewportCount = 1;
        viewportState.pViewports = &viewport;
        viewportState.scissorCount = 1;
        viewportState.pScissors = &scissor;

        VkPipelineRasterizationStateCreateInfo rasterizer{};
        rasterizer.sType = VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO;
        rasterizer.depthClampEnable = VK_FALSE;
        rasterizer.rasterizerDiscardEnable = VK_FALSE;
        rasterizer.polygonMode = VK_POLYGON_MODE_FILL;
        rasterizer.lineWidth = 1.0f;
        rasterizer.cullMode = VK_CULL_MODE_BACK_BIT;
        rasterizer.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;
        rasterizer.depthBiasEnable = VK_FALSE;

        VkPipelineMultisampleStateCreateInfo multisampling{};
        multisampling.sType = VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO;
        multisampling.sampleShadingEnable = VK_FALSE;
        multisampling.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

        VkPipelineDepthStencilStateCreateInfo depthStencil{};
        depthStencil.sType = VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO;
        depthStencil.depthTestEnable = VK_TRUE;
        depthStencil.depthWriteEnable = VK_TRUE;
        depthStencil.depthCompareOp = VK_COMPARE_OP_LESS;
        depthStencil.depthBoundsTestEnable = VK_FALSE;
        depthStencil.stencilTestEnable = VK_FALSE;

        VkPipelineColorBlendAttachmentState colorBlendAttachment{};
        colorBlendAttachment.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
        colorBlendAttachment.blendEnable = VK_FALSE;

        VkPipelineColorBlendStateCreateInfo colorBlending{};
        colorBlending.sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO;
        colorBlending.logicOpEnable = VK_FALSE;
        colorBlending.logicOp = VK_LOGIC_OP_COPY;
        colorBlending.attachmentCount = 1;
        colorBlending.pAttachments = &colorBlendAttachment;
        colorBlending.blendConstants[0] = 0.0f;
        colorBlending.blendConstants[1] = 0.0f;
        colorBlending.blendConstants[2] = 0.0f;
        colorBlending.blendConstants[3] = 0.0f;

        VkGraphicsPipelineCreateInfo pipelineInfo{};
        pipelineInfo.sType = VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO;
        pipelineInfo.stageCount = 2;
        pipelineInfo.pStages = shaderStages;
        pipelineInfo.pVertexInputState = &vertexInputInfo;
        pipelineInfo.pInputAssemblyState = &inputAssembly;
        pipelineInfo.pViewportState = &viewportState;
        pipelineInfo.pRasterizationState = &rasterizer;
        pipelineInfo.pMultisampleState = &multisampling;
        pipelineInfo.pDepthStencilState = &depthStencil;
        pipelineInfo.pColorBlendState = &colorBlending;
        pipelineInfo.layout = layout;
        pipelineInfo.renderPass = render_pass;
        pipelineInfo.subpass = 0;
        pipelineInfo.basePipelineHandle = VK_NULL_HANDLE;

        NVVK_CHECK(vkCreateGraphicsPipelines(device, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &pipeline));

        vkDestroyShaderModule(device, fragShaderModule, nullptr);
        vkDestroyShaderModule(device, vertShaderModule, nullptr);
    }

    MeshPipeline(MeshPipeline&&) = delete;

    ~MeshPipeline()
    {
        vkDestroyPipeline(device, pipeline, nullptr);
        vkDestroyPipelineLayout(device, layout, nullptr);
    }

    operator VkPipeline() const
    {
        return pipeline;
    }
};



} // end anonymous namespace.



namespace myricube {

struct RendererVk :
    RendererLogic<MeshEntry, MeshStaging, RaycastEntry, RaycastStaging>
{
    // The order of these is very important.
    ScopedContext ctx;
    ScopedSurface glfw_surface;
    FrameManager frame_manager;
    ScopedRenderPass render_pass;
    std::unique_ptr<MeshPipeline> p_mesh_pipeline;

    // Set in begin_frame.
    int width, height;

    VkCommandBuffer frame_cmd_buffer = VK_NULL_HANDLE;
    nvvk::SwapChainImage current_swap_image;

    RendererVk(RenderThread* thread, RenderArgs args) :
        RendererLogic<MeshEntry, MeshStaging, RaycastEntry, RaycastStaging>(
            thread,
            args),
        ctx(),
        glfw_surface(
            ctx.m_instance,
            args.p_window->get_glfw_window()),
        frame_manager(
            ctx,
            glfw_surface.surface,
            glfw_surface.initial_width,
            glfw_surface.initial_height),
        render_pass(
            ctx.m_device)
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
        if (p_mesh_pipeline == nullptr
            or p_mesh_pipeline->width != width
            or p_mesh_pipeline->height != height)
        {
            p_mesh_pipeline.reset(new MeshPipeline(
                ctx.m_device,
                render_pass,
                width, height));
        }
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
