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
#include <mutex>
#include <stdexcept>
#include <string.h>
#include <vulkan/vulkan.h>

#include "nvvk/context_vk.hpp"
#include "FrameManager.hpp"

using namespace myricube;
using namespace akeley;

namespace myricube { struct RendererVk; }

namespace {

// Set by RendererVk constructor so other constructors/destructors can
// access it. (But note this only occurs after RendererVk is fully
// constructed).
thread_local RendererVk* renderer = nullptr;

// RendererVk is responsible for constructing MeshEntry and so on.
// Need these helpers for now as RendererVk is only forward-declared.
struct MeshEntry;
void constructor(MeshEntry*);
void destructor(MeshEntry*);

// These should be widely supported, but I can detect support if
// needed later.
constexpr auto swap_chain_image_format = VK_FORMAT_B8G8R8A8_UNORM;
constexpr auto depth_format = VK_FORMAT_D24_UNORM_S8_UINT;

// Handle for GPU resources for the mesh of one chunk group.
struct MeshEntry
{
    // Memory and buffer for the voxel instance list.
    // TODO: Share memory among multiple MeshEntry.
    VkDeviceMemory memory;
    VkBuffer buffer;

    // Persistent coherent mapping of buffer.
    MappedGroupMesh* map;

    // Data needed to draw chunk[z][y][x] within this chunk group.
    ChunkDrawData draw_data[edge_chunks][edge_chunks][edge_chunks];

    // Size in bytes of the vbo's data store on the GPU.
    static constexpr GLsizeiptr vbo_bytes = sizeof(MappedGroupMesh);

    // Return the offset (in number of MeshVoxelVertex's, not bytes)
    // into the VBO where the data from mesh_array[z][y][x] is copied
    // into.
    static unsigned vert_offset(unsigned x, unsigned y, unsigned z)
    {
        assert(x < edge_chunks and y < edge_chunks and z < edge_chunks);
        unsigned chunk_idx = x + y*edge_chunks + z*edge_chunks*edge_chunks;
        return chunk_idx * chunk_max_verts;
    }

    // Same as above, but return offset as count of bytes.
    static VkDeviceSize byte_offset(unsigned x, unsigned y, unsigned z)
    {
        auto vert_sz = VkDeviceSize(sizeof(MeshVoxelVertex));
        VkDeviceSize off = vert_offset(x, y, z) * vert_sz;
        assert(size_t(off) < vbo_bytes);
        return off;
    }

    MeshEntry() { constructor(this); }
    MeshEntry(MeshEntry&&) = delete;
    ~MeshEntry() { destructor(this); }
};

void swap(MeshEntry& left, MeshEntry& right) noexcept
{
    using std::swap;
    swap(left.memory, right.memory);
    swap(left.buffer, right.buffer);
    swap(left.map, right.map);
    swap(left.draw_data, right.draw_data);
}

struct MeshStaging
{
    MeshEntry entry;
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



// camelCase functions copied with some modifications from
// vulkan-tutorial.com; they're mostly exempt from my 80-column limit
// for now.
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
    NVVK_CHECK(vkCreateShaderModule(device, &createInfo, nullptr, &shaderModule));

    return shaderModule;
}

uint32_t findMemoryType(
    VkPhysicalDevice physicalDevice,
    uint32_t typeFilter,
    VkMemoryPropertyFlags properties)
{
    VkPhysicalDeviceMemoryProperties memProperties;
    vkGetPhysicalDeviceMemoryProperties(physicalDevice, &memProperties);

    for (uint32_t i = 0; i < memProperties.memoryTypeCount; i++) {
        if ((typeFilter & (1 << i)) && (memProperties.memoryTypes[i].propertyFlags & properties) == properties) {
            return i;
        }
    }

    throw std::runtime_error("failed to find suitable memory type!");
}

// Create a buffer with dedicated device memory.
void createBuffer(
    VkPhysicalDevice physicalDevice,
    VkDevice device,
    VkDeviceSize size,
    VkBufferUsageFlags usage,
    VkMemoryPropertyFlags properties,
    VkBuffer& buffer,
    VkDeviceMemory& bufferMemory)
{
    VkBufferCreateInfo bufferInfo{};
    bufferInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
    bufferInfo.size = size;
    bufferInfo.usage = usage;
    bufferInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

    NVVK_CHECK(vkCreateBuffer(device, &bufferInfo, nullptr, &buffer));

    VkMemoryRequirements memRequirements;
    vkGetBufferMemoryRequirements(device, buffer, &memRequirements);

    VkMemoryAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocInfo.allocationSize = memRequirements.size;
    allocInfo.memoryTypeIndex = findMemoryType(physicalDevice, memRequirements.memoryTypeBits, properties);

    NVVK_CHECK(vkAllocateMemory(device, &allocInfo, nullptr, &bufferMemory));

    vkBindBufferMemory(device, buffer, bufferMemory, 0);
}

// Create an image with dedicated device memory.
void createImage(
    VkPhysicalDevice physicalDevice,
    VkDevice device,
    uint32_t width, uint32_t height,
    VkFormat format,
    VkImageTiling tiling,
    VkImageUsageFlags usage,
    VkMemoryPropertyFlags properties,
    VkImage& image,
    VkDeviceMemory& imageMemory)
{
    VkImageCreateInfo imageInfo{};
    imageInfo.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
    imageInfo.imageType = VK_IMAGE_TYPE_2D;
    imageInfo.extent.width = width;
    imageInfo.extent.height = height;
    imageInfo.extent.depth = 1;
    imageInfo.mipLevels = 1;
    imageInfo.arrayLayers = 1;
    imageInfo.format = format;
    imageInfo.tiling = tiling;
    imageInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    imageInfo.usage = usage;
    imageInfo.samples = VK_SAMPLE_COUNT_1_BIT;
    imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

    if (vkCreateImage(device, &imageInfo, nullptr, &image) != VK_SUCCESS) {
        throw std::runtime_error("failed to create image!");
    }

    VkMemoryRequirements memRequirements;
    vkGetImageMemoryRequirements(device, image, &memRequirements);

    VkMemoryAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocInfo.allocationSize = memRequirements.size;
    allocInfo.memoryTypeIndex = findMemoryType(physicalDevice, memRequirements.memoryTypeBits, properties);

    if (vkAllocateMemory(device, &allocInfo, nullptr, &imageMemory) != VK_SUCCESS) {
        throw std::runtime_error("failed to allocate image memory!");
    }

    vkBindImageMemory(device, image, imageMemory, 0);
}

VkImageView createImageView(
    VkDevice device,
    VkImage image,
    VkFormat format,
    VkImageAspectFlags aspectFlags)
{
    VkImageViewCreateInfo viewInfo{};
    viewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
    viewInfo.image = image;
    viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
    viewInfo.format = format;
    viewInfo.subresourceRange.aspectMask = aspectFlags;
    viewInfo.subresourceRange.baseMipLevel = 0;
    viewInfo.subresourceRange.levelCount = 1;
    viewInfo.subresourceRange.baseArrayLayer = 0;
    viewInfo.subresourceRange.layerCount = 1;

    VkImageView imageView;
    NVVK_CHECK(vkCreateImageView(device, &viewInfo, nullptr, &imageView));
    return imageView;
}



// Manager for framebuffers + shared depth buffer, one per swap chain image.
struct Framebuffers
{
    // From nvvk::SwapChain::getChangeID().  Basically, if this
    // doesn't match that of nvvk::SwapChain, the swap chain has been
    // re-created, and we need to re-create the framebuffers here to
    // match.
    uint32_t last_change_id;

    // Borrowed device pointers and render pass.
    VkPhysicalDevice physical_device;
    VkDevice device;
    VkRenderPass render_pass;

    // Shared depth buffer and its memory and ImageView.
    VkDeviceMemory depth_memory;
    VkImage depth_image;
    VkImageView depth_view;

    // framebuffer[i] is the framebuffer for swap image i, as you'd
    // expect.  This is cleared to indicate when this class is in an
    // unitinialized state.
    std::vector<VkFramebuffer> framebuffers;
    bool initialized() const { return !framebuffers.empty(); }

    Framebuffers(
        VkPhysicalDevice physical_device_,
        VkDevice dev_,
        VkRenderPass render_pass_) :
    physical_device(physical_device_),
    device(dev_),
    render_pass(render_pass_) { }

    ~Framebuffers() { destructor(); }

    // Check the swap chain and recreate now if needed (now = no
    // synchronization done; note however that we can rely on
    // FrameManager to wait on the main thread queue to idle before
    // re-creating a swap chain).
    void recreate_now_if_needed(nvvk::SwapChain& swap_chain) noexcept
    {
        if (initialized() and swap_chain.getChangeID() == last_change_id) {
            return;
        }

        // Destroy old resources.
        destructor();

        // Make depth buffer.
        createImage(
            physical_device,
            device,
            swap_chain.getUpdateWidth(),
            swap_chain.getUpdateHeight(),
            depth_format,
            VK_IMAGE_TILING_OPTIMAL,
            VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
            depth_image,
            depth_memory);

        depth_view = createImageView(
            device, depth_image, depth_format, VK_IMAGE_ASPECT_DEPTH_BIT);

        // Make a framebuffer for every swap chain image.
        size_t image_count = swap_chain.getImageCount();
        framebuffers.resize(image_count);
        for (size_t i = 0; i < image_count; ++i) {
            std::array<VkImageView, 2> attachments = {
                swap_chain.getImageView(i),
                depth_view
            };

            VkFramebufferCreateInfo framebufferInfo{};
            framebufferInfo.sType = VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO;
            framebufferInfo.renderPass = render_pass;
            framebufferInfo.attachmentCount =
                static_cast<uint32_t>(attachments.size());
            framebufferInfo.pAttachments = attachments.data();
            framebufferInfo.width = swap_chain.getUpdateWidth();
            framebufferInfo.height = swap_chain.getUpdateHeight();
            framebufferInfo.layers = 1;

            NVVK_CHECK(vkCreateFramebuffer(
                device, &framebufferInfo, nullptr, &framebuffers.at(i)));
        }

        last_change_id = swap_chain.getChangeID();
    }

    void destructor()
    {
        if (initialized()) {
            vkDestroyImageView(device, depth_view, nullptr);
            vkDestroyImage(device, depth_image, nullptr);
            vkFreeMemory(device, depth_memory, nullptr);
            for (VkFramebuffer fb : framebuffers) {
                vkDestroyFramebuffer(device, fb, nullptr);
            }
            framebuffers.clear();
        }
        assert(!initialized());
    }

    VkFramebuffer operator[] (size_t i) const
    {
        return framebuffers.at(i);
    }
};



struct PushConstant
{
    glm::mat4 mvp;
    glm::vec4 eye_relative_group_origin; // w unused.
    int32_t flags;
    int32_t far_plane_squared;
    int32_t raycast_thresh_squared;
};
static_assert(sizeof(PushConstant) <= 128);



// Pipeline for mesh rendering. All we really have to do here is load
// the vert/frag shader and hook up the instanced voxels
// array. Everything else is just endless boilerplate.
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
        pushConstantRange.stageFlags = VK_SHADER_STAGE_ALL_GRAPHICS;
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
    VkDevice dev;
    ScopedSurface glfw_surface;
    FrameManager frame_manager;
    ScopedRenderPass render_pass;
    Framebuffers framebuffers;
    std::unique_ptr<MeshPipeline> p_mesh_pipeline;

    // Set in begin_frame.
    int width, height;

    // Graphics+Present queue used by the main thread, and its family.
    VkQueue main_queue;
    uint32_t main_queue_family;

    // Transfer-only queue shared by the worker threads, plus a mutex
    // synchronizing it.
    VkQueue worker_queue;
    uint32_t worker_queue_family;
    std::mutex worker_queue_mutex;

    // From FrameManager: primary command buffer that will be
    // submitted at frame-end, and the current swap chain image. Both
    // owned by the main thread.
    VkCommandBuffer frame_cmd_buffer = VK_NULL_HANDLE;
    nvvk::SwapChainImage current_swap_image;

    RendererVk(RenderThread* thread, RenderArgs args) :
        RendererLogic<MeshEntry, MeshStaging, RaycastEntry, RaycastStaging>(
            thread,
            args),
        ctx(),
        dev(ctx.m_device),
        glfw_surface(
            ctx.m_instance,
            args.p_window->get_glfw_window()),
        frame_manager(
            ctx,
            glfw_surface.surface,
            glfw_surface.initial_width,
            glfw_surface.initial_height),
        render_pass(
            dev),
        framebuffers(
            ctx.m_physicalDevice,
            dev,
            render_pass)
    {
        // As promised.
        renderer = this;

        // Steal the queues from the ever-helpful nvvk::Context.
        // ctx.m_queueGCT is the same queue FrameManager is using.
        main_queue =        ctx.m_queueGCT;
        main_queue_family = ctx.m_queueGCT;
        worker_queue =        ctx.m_queueT;
        worker_queue_family = ctx.m_queueT;

        if (main_queue == VK_NULL_HANDLE) {
            throw std::runtime_error("Failed to find graphics+present queue");
        }
        if (worker_queue == VK_NULL_HANDLE) {
            throw std::runtime_error("Failed to find transfer-only queue");
        }
    }

    void begin_frame() override
    {
        // Begin the frame, recreating the swap chain (hidden in
        // beginFrame) and framebuffers if needed.
        glfwGetFramebufferSize(glfw_surface.p_window, &width, &height);
        frame_manager.beginFrame(
            &frame_cmd_buffer,
            &current_swap_image,
            width, height);
        framebuffers.recreate_now_if_needed(frame_manager.getSwapChain());

        VkClearColorValue clear_color;
        clear_color.float32[0] = 0.0f;
        clear_color.float32[1] = 0.0f;
        clear_color.float32[2] = 0.0f;
        clear_color.float32[3] = 1.0f;
        std::array<VkClearValue, 2> clear_values{};
        clear_values[0].color = clear_color;            // Color attachment
        clear_values[1].depthStencil = { 1.0f, 0 };     // Depth attachment

        // Start the render pass, drawing into the correct framebuffer
        // for the current swap chain image.
        VkRenderPassBeginInfo begin_info {
            VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO,
            nullptr,
            render_pass,
            framebuffers[current_swap_image.index],
            { { 0, 0 }, { uint32_t(width), uint32_t(height) } },
            clear_values.size(),
            clear_values.data()
        };
        vkCmdBeginRenderPass(
            frame_cmd_buffer, &begin_info, VK_SUBPASS_CONTENTS_INLINE);
    }


    /* MESH RENDERING IMPLEMENTATION */

    // Add to the primary command buffer commands to draw the list
    // of chunk groups with the instanced-voxel rendering method.
    // At this point, we are in the render pass begin_frame started.
    void draw_mesh_entries(
        const std::vector<std::pair<MeshEntry*, glm::ivec3>>& entries) override
    {
        // Since I don't know how to use dynamic state, we have to recreate
        // the pipeline if the framebuffer size changed.
        if (p_mesh_pipeline == nullptr
            or p_mesh_pipeline->width != width
            or p_mesh_pipeline->height != height)
        {
            p_mesh_pipeline.reset(new MeshPipeline(
                dev,
                render_pass,
                width, height));
        }

        PushConstant push_constant{};
        glm::mat4 vp = transforms.residue_vp_matrix;
        glm::vec3 eye_residue = transforms.eye_residue;
        glm::ivec3 eye_group = transforms.eye_group;

        // Bind the mesh-drawing pipeline.
        vkCmdBindPipeline(
            frame_cmd_buffer,
            VK_PIPELINE_BIND_POINT_GRAPHICS,
            *p_mesh_pipeline);

        for (auto pair : entries) {
            MeshEntry& entry = *pair.first;
            glm::ivec3 group_coord = pair.second;

            // Push the push constant. We need to update the MVP.
            // The view matrix only takes into account the eye's residue
            // coordinates, so the model position of the group actually
            // needs to be shifted by the eye's group coord.
            glm::vec3 model_offset = glm::vec3(group_coord - eye_group)
                                   * float(group_size);
            glm::mat4 m = glm::translate(glm::mat4(1.0f), model_offset);
            push_constant.mvp = vp * m;
            vkCmdPushConstants(
                frame_cmd_buffer,
                p_mesh_pipeline->layout,
                VK_SHADER_STAGE_ALL_GRAPHICS,
                0, sizeof(PushConstant), &push_constant);

            // Bind the instanced voxels buffer for this group in
            // binding point 0.
            VkDeviceSize offsets{0};
            vkCmdBindVertexBuffers(
                frame_cmd_buffer, 0, 1, &entry.buffer, &offsets);
            
            // Draw every chunk. Each has its voxels starting at a
            // different index in the buffer.
            for (int z = 0; z < edge_chunks; ++z) {
            for (int y = 0; y < edge_chunks; ++y) {
            for (int x = 0; x < edge_chunks; ++x) {
                auto draw_data = entry.draw_data[z][y][x];
                PackedAABB aabb = draw_data.aabb;
                if (decide_chunk(group_coord, aabb) != draw_mesh) {
                    continue;
                }

                if (draw_data.vert_count == 0) continue;
                // Draw with instance rendering (36 verts per voxel).
                uint32_t instances = uint32_t(draw_data.vert_count);
                uint32_t base = MeshEntry::vert_offset(x, y, z);
                vkCmdDraw(frame_cmd_buffer, 36, instances, 0, base);
            }
            }
            }
        }
    }
    
    // Construct the given MeshEntry. Only called by main thread. TODO
    void constructor(MeshEntry* entry)
    {
        VkDeviceSize bytes = sizeof(MappedGroupMesh);
        createBuffer(
            ctx.m_physicalDevice,
            dev,
            bytes,
            VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
            VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
            entry->buffer,
            entry->memory);
        void* p_void;
        vkMapMemory(dev, entry->memory, 0, bytes, 0, &p_void);
        entry->map = static_cast<MappedGroupMesh*>(p_void);
    }
    
    // Release resources of MeshEntry.  Only called by main thread. TODO
    void destructor(MeshEntry* entry)
    {
        vkUnmapMemory(dev, entry->memory);
        vkDestroyBuffer(dev, entry->buffer, nullptr);
        vkFreeMemory(dev, entry->memory, nullptr);
    }

    // Convert chunk group's chunks to meshes and pack into the
    // memory-mapped array of chunks.
    void worker_stage(
        MeshStaging* staging, const BinChunkGroup* group_ptr) override
    {
        MeshEntry* entry = &staging->entry;
        for (int zL = 0; zL < edge_chunks; ++zL) {
        for (int yL = 0; yL < edge_chunks; ++yL) {
        for (int xL = 0; xL < edge_chunks; ++xL) {
            fill_chunk_mesh(&entry->map->chunks[zL][yL][xL],
                            &entry->draw_data[zL][yL][xL],
                            group_ptr->chunk_array[zL][yL][xL],
                            glm::ivec3(xL, yL, zL) * chunk_size);
        }
        }
        }
    }

    bool swap_in(
        MeshStaging* staging,
        std::unique_ptr<MeshEntry>* p_uptr_entry) override
    {
        if (*p_uptr_entry == nullptr) p_uptr_entry->reset(new MeshEntry());
        swap(staging->entry, **p_uptr_entry);
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
        // End the one render pass there is.
        vkCmdEndRenderPass(frame_cmd_buffer);

        // Ends cmd buffer recording, submits it, and presents the
        // rendered image.
        frame_manager.endFrame(frame_cmd_buffer);
    }

    void wait_idle() override
    {
        vkDeviceWaitIdle(dev);
    }
};

std::shared_ptr<RendererBase> RendererVk_Factory(
    RenderThread* thread,
    RenderArgs args)
{
    return std::shared_ptr<RendererBase>(new RendererVk(thread, args));
}

} // end namespace myricube

namespace {

void constructor(MeshEntry* entry) { renderer->constructor(entry); }
void destructor(MeshEntry* entry) { renderer->destructor(entry); }

}