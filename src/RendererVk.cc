// Vulkan implementation of RendererLogic.

#define GLFW_INCLUDE_NONE
#define GLFW_INCLUDE_VULKAN

#include "myricube.hh"
#include "RendererLogic.hh"

#include "GLFW/glfw3.h"

#include <cassert>
#include <errno.h>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <string.h>
#include <vulkan/vulkan.h>

#include "nvvk/context_vk.hpp"
#include "EnvVar.hh"
#include "FrameManager.hpp"

#include "../myricube-data/vk/PushConstant.glsl"

using namespace myricube;
using namespace vkdo;

namespace myricube {
    struct RendererVk;
    bool vulkan_compiled_in = true;
}

namespace {

EnvVar64 validation_enabled("myricube_validation", 0);

// Set by RendererVk constructor so other constructors/destructors can
// access it. (But note this only occurs after RendererVk is fully
// constructed).
thread_local RendererVk* thread_local_renderer = nullptr;

constexpr float min_depth = 0.0f, max_depth = 1.0f;

// RendererVk is responsible for constructing MeshEntry and so on.
// Need these helpers for now as RendererVk is only forward-declared.
struct MeshEntry;
void constructor(MeshEntry*);
void destructor(MeshEntry*);
struct RaycastEntry;
void constructor(RaycastEntry*);
void destructor(RaycastEntry*);
struct RaycastStaging;
void constructor(RaycastStaging*);
void destructor(RaycastStaging*);

// These should be widely supported, but I can detect support if
// needed later. No stencil buffer for now.
constexpr auto swap_chain_image_format = VK_FORMAT_B8G8R8A8_UNORM;
constexpr auto depth_format = VK_FORMAT_D32_SFLOAT;

// Image format for 3D images used to store chunk groups' voxels.
// Need to check that it matches the bit assignments.
constexpr auto chunk_group_voxels_image_format = VK_FORMAT_R8G8B8A8_UNORM;
static_assert(red_shift == 0);
static_assert(green_shift == 8);
static_assert(blue_shift == 16);
// Assume little endian for now, the endian wars have finally ended.

constexpr VkImageSubresourceRange color_range {
    VK_IMAGE_ASPECT_COLOR_BIT,
    0, VK_REMAINING_MIP_LEVELS,
    0, VK_REMAINING_ARRAY_LAYERS };

// Handle for GPU resources for the mesh of one chunk group.
struct MeshEntry
{
    // Buffer for the voxel instance list.
    VkBuffer buffer = VK_NULL_HANDLE;

    // Persistent coherent mapping of buffer.
    MappedGroupMesh* map = nullptr;

    // Data needed to draw chunk[z][y][x] within this chunk group.
    ChunkDrawData draw_data[edge_chunks][edge_chunks][edge_chunks];

    // Make sure you fix the swap function if you add more stuff.

    // Size in bytes of the vbo's data store on the GPU.
    static constexpr size_t vbo_bytes = sizeof(MappedGroupMesh);

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

    MeshEntry(VkBuffer buffer_, MappedGroupMesh* map_)
    {
        buffer = buffer_;
        map = map_;
    }

    MeshEntry(MeshEntry&& other) = delete;

    ~MeshEntry() { destructor(this); }
};

void swap(MeshEntry& left, MeshEntry& right) noexcept
{
    using std::swap;
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
    // VBO storing packed AABBs, and a memory map of it.
    VkBuffer aabb_buffer = VK_NULL_HANDLE;
    AABBs* aabb_map = nullptr;

    // 3D image storing voxels of this chunk group.
    VkImage voxels_image = VK_NULL_HANDLE;
    VkImageView voxels_view = VK_NULL_HANDLE;

    // Descriptor set referring to above 3D image.
    VkDescriptorSet voxels_descriptor = VK_NULL_HANDLE;

    RaycastEntry() { constructor(this); }

    RaycastEntry(RaycastEntry&&) = delete;

    ~RaycastEntry() { destructor(this); }
};



struct MappedChunks
{
    BinChunk chunks[edge_chunks][edge_chunks][edge_chunks];
};

// Staging buffer for the async cache of raycast chunk groups.
struct RaycastStaging
{
    // RaycastEntry to be created.
    std::unique_ptr<RaycastEntry> entry;

    // Buffer storing the staging buffer, and a memory map of it.
    // The chunks are stored in [z][y][x] order as usual.
    VkBuffer chunk_buffer = VK_NULL_HANDLE;
    MappedChunks* chunk_map = nullptr;

    // Fence for the main thread to wait on to indicate the staging
    // buffer -> image copy is done, and we can swap the RaycastEntry
    // into the main cache.
    std::shared_ptr<VkFence> p_fence = nullptr;

    // Used in swap_in RaycastStaging.
    bool command_recorded = false;

    RaycastStaging() { constructor(this); }

    RaycastStaging(RaycastStaging&&) = delete;

    ~RaycastStaging() { destructor(this); }
};

// nvvk::Context with proper constructor and destructor.
struct ScopedContext : nvvk::Context
{
    ScopedContext()
    {
        // Add the swapchain extensions.
        nvvk::ContextCreateInfo context_info(validation_enabled);
        context_info.verboseCompatibleDevices = false;
        context_info.verboseUsed = true;
        context_info.verboseAvailable = false;
        context_info.apiMajor = 1;
        context_info.apiMinor = 1;
        uint32_t extension_count;
        const char** extensions =
            glfwGetRequiredInstanceExtensions(&extension_count);
        if (extensions == nullptr) {
            throw std::runtime_error("Could not get glfw vk extensions.");
        }
        for (uint32_t i = 0; i < extension_count; ++i) {
            context_info.addInstanceExtension(extensions[i]);
        }
        context_info.addDeviceExtension(VK_KHR_SWAPCHAIN_EXTENSION_NAME);

        // Initialize the nvvk context.
        init(context_info);

        // Warn if validation disabled.
        if (!validation_enabled) {
            fprintf(stderr,
            #if !defined(MYRICUBE_WINDOWS)
                "\x1b[35m\x1b[1m"
            #endif
                "Note for developers: Vulkan Validation Layers disabled!\n"
            #if !defined(MYRICUBE_WINDOWS)
                "\x1b[0m"
            #endif
                "(enable with environment variable myricube_validation=1)\n"
            );
        }
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
std::vector<char> readFile(const filename_string& filename)
{
#ifdef MYRICUBE_WINDOWS
    std::string ascii_filename;
    ascii_filename.reserve(filename.size());
    for (filename_char c : filename) {
        if (c <= 127) ascii_filename.push_back(char(c));
        else throw std::runtime_error(
            "TODO: support non-ASCII Vulkan shader path"
            " (call me and complain)");
        // I should just switch to _wfopen.
    }
#else
    const std::string& ascii_filename = filename;
#endif

    std::ifstream file(ascii_filename, std::ios::ate | std::ios::binary);

    if (!file.is_open()) {
        const char* err = strerror(errno);
        #ifdef MYRICUBE_WINDOWS
            fprintf(stderr,
                "Failed to open file: %ls (%s)\n", filename.c_str(), err);
            throw std::runtime_error("Failed to open file");
        #else
            throw std::runtime_error("failed to open file: " + filename
                + " (" + err + ")");
        #endif
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

template <typename T>
struct MappedBuffer
{
    VkBuffer buffer;
    T* map;
};

// Create buffer_count many buffers of buffer_bytes bytes each that
// all share one VK memory allocation, which is returned through
// *p_buffer_memory. The memory is mapped to the host, and each
// buffer is paired with a pointer to its portion of the mapping.
template <typename T>
std::vector<MappedBuffer<T>> create_mapped_buffer_array(
    VkPhysicalDevice physical_device,
    VkDevice device,
    VkDeviceSize buffer_bytes,
    size_t buffer_count,
    VkBufferUsageFlags usage,
    VkMemoryPropertyFlags properties,
    VkDeviceMemory* p_buffer_memory) noexcept
{
    assert(buffer_count != 0);

    VkBufferCreateInfo bufferInfo{};
    bufferInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
    bufferInfo.size = buffer_bytes;
    bufferInfo.usage = usage;
    bufferInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

    std::vector<MappedBuffer<T>> array(buffer_count);
    for (auto& map : array) {
        NVVK_CHECK(vkCreateBuffer(device, &bufferInfo, nullptr, &map.buffer));
    }

    VkMemoryRequirements requirements;
    vkGetBufferMemoryRequirements(device, array[0].buffer, &requirements);
    VkDeviceSize stride = (requirements.size + requirements.alignment - 1)
                        & ~(requirements.alignment - 1);

    VkMemoryAllocateInfo allocInfo{};
    VkDeviceSize alloc_size = stride * buffer_count;
    allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocInfo.allocationSize = alloc_size;
    allocInfo.memoryTypeIndex = findMemoryType(
        physical_device, requirements.memoryTypeBits, properties);

    NVVK_CHECK(vkAllocateMemory(device, &allocInfo, nullptr, p_buffer_memory));
    VkDeviceMemory memory = *p_buffer_memory;

    void* p_void;
    NVVK_CHECK(vkMapMemory(device, memory, 0, alloc_size, 0, &p_void));
    char* memory_mapping = static_cast<char*>(p_void);

    for (size_t i = 0; i < buffer_count; ++i) {
        VkDeviceSize off = stride * i;
        NVVK_CHECK(vkBindBufferMemory(device, array[i].buffer, memory, off));
        array[i].map = reinterpret_cast<T*>(memory_mapping + off);
    }
    return array;
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

    NVVK_CHECK(vkCreateImage(device, &imageInfo, nullptr, &image));

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



// Create a vector of 3D images suitable for storing a chunk group's voxels,
// paired with a corresponding image view. All voxels share one memory
// allocation, returned through *p_memory.
std::vector<std::pair<VkImage, VkImageView>> create_chunk_group_images(
    VkPhysicalDevice physicalDevice,
    VkDevice device,
    size_t image_count,
    VkFormat format,
    VkImageTiling tiling,
    VkImageUsageFlags usage,
    VkMemoryPropertyFlags properties,
    VkDeviceMemory* p_memory)
{
    VkImageCreateInfo imageInfo{};
    imageInfo.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
    imageInfo.imageType = VK_IMAGE_TYPE_3D;
    imageInfo.extent.width = group_size;
    imageInfo.extent.height = group_size;
    imageInfo.extent.depth = group_size;
    imageInfo.mipLevels = 1;
    imageInfo.arrayLayers = 1;
    imageInfo.format = format;
    imageInfo.tiling = tiling;
    imageInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    imageInfo.usage = usage;
    imageInfo.samples = VK_SAMPLE_COUNT_1_BIT;
    imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

    assert(image_count != 0);
    std::vector<std::pair<VkImage, VkImageView>> result(image_count);

    VkImage* p_image;
    for (size_t i = 0; i < image_count; ++i) {
        p_image = &result[i].first;
        NVVK_CHECK(vkCreateImage(device, &imageInfo, nullptr, p_image));
    }

    VkMemoryRequirements requirements;
    vkGetImageMemoryRequirements(device, *p_image, &requirements);

    VkDeviceSize stride = (requirements.size + requirements.alignment - 1)
                        & ~(requirements.alignment - 1);

    VkMemoryAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocInfo.allocationSize = stride * VkDeviceSize(image_count);
    allocInfo.memoryTypeIndex = findMemoryType(
        physicalDevice, requirements.memoryTypeBits, properties);

    NVVK_CHECK(vkAllocateMemory(device, &allocInfo, nullptr, p_memory));

    VkImageViewCreateInfo view_info{};
    view_info.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
    view_info.image = VK_NULL_HANDLE;
    view_info.viewType = VK_IMAGE_VIEW_TYPE_3D;
    view_info.format = format;
    view_info.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    view_info.subresourceRange.baseMipLevel = 0;
    view_info.subresourceRange.levelCount = 1;
    view_info.subresourceRange.baseArrayLayer = 0;
    view_info.subresourceRange.layerCount = 1;

    for (VkDeviceSize i = 0; i < VkDeviceSize(image_count); ++i) {
        VkImage image = result[i].first;
        NVVK_CHECK(vkBindImageMemory(device, image, *p_memory, stride * i));
        view_info.image = image;
        NVVK_CHECK(vkCreateImageView(
            device, &view_info, nullptr, &result[i].second));
    }

    return result;
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
            swap_chain.getWidth(),
            swap_chain.getHeight(),
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
            framebufferInfo.width = swap_chain.getWidth();
            framebufferInfo.height = swap_chain.getHeight();
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



// Dynamic viewport/scissor pipeline for mesh rendering. All we really
// have to do here is load the vert/frag shader and hook up the
// instanced voxels array. Everything else is just endless
// boilerplate.
struct MeshPipeline
{
    // We manage these.
    VkPipeline pipeline;
    VkPipelineLayout layout;

    // Borrowed pointer.
    VkDevice device;

    MeshPipeline(VkDevice dev_, VkRenderPass render_pass)
    {
        device = dev_;

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
        viewport.width = 1.0f;
        viewport.height = 1.0f;
        viewport.minDepth = min_depth;
        viewport.maxDepth = max_depth;

        VkRect2D scissor{};
        scissor.offset = {0, 0};
        scissor.extent = VkExtent2D{ 1, 1 };

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
        pipelineInfo.flags = VK_DYNAMIC_STATE_VIEWPORT | VK_DYNAMIC_STATE_SCISSOR;
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



// Raycast pipeline. We have to draw AABBs with instanced rendering
// and set up a pipeline layout for feeding a push constant plus a 3D
// voxels image.
struct RaycastPipeline
{
    // We manage these.
    VkPipeline pipeline;
    VkPipelineLayout pipeline_layout;
    VkDescriptorSetLayout set_layout;

    // Borrowed pointer.
    VkDevice device;

    static constexpr uint32_t pool_sets = 1024;
    static constexpr VkDescriptorPoolSize pool_size {
        VK_DESCRIPTOR_TYPE_STORAGE_IMAGE,
        pool_sets };

    // A VkDescriptorPoolCreateInfo that describes creating a pool
    // with enough room for pool_sets-many descriptors for the raycast
    // pipeline.
    static constexpr VkDescriptorPoolCreateInfo descriptor_pool_info {
        VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO,
        nullptr,
        0,
        pool_sets,
        1,
        &pool_size };

    RaycastPipeline(VkDevice dev_, VkRenderPass render_pass)
    {
        device = dev_;

        // Boilerplate from vulkan-tutorial.com with some modifications.

        // Only binding is the 3D image holding voxel data.
        VkDescriptorSetLayoutBinding voxels_binding {
            0,
            VK_DESCRIPTOR_TYPE_STORAGE_IMAGE,
            1,
            VK_SHADER_STAGE_FRAGMENT_BIT,
            nullptr };

        VkDescriptorSetLayoutCreateInfo set_info {
            VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO,
            nullptr,
            0,
            1,
            &voxels_binding };

        NVVK_CHECK(vkCreateDescriptorSetLayout(
            device, &set_info, nullptr, &set_layout));

        VkPushConstantRange pushConstantRange{};
        pushConstantRange.stageFlags = VK_SHADER_STAGE_ALL_GRAPHICS;
        pushConstantRange.offset = 0;
        pushConstantRange.size = sizeof(PushConstant);

        VkPipelineLayoutCreateInfo pipelineLayoutInfo{};
        pipelineLayoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        pipelineLayoutInfo.setLayoutCount = 1;
        pipelineLayoutInfo.pSetLayouts = &set_layout;
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushConstantRange;

        NVVK_CHECK(vkCreatePipelineLayout(
            device, &pipelineLayoutInfo, nullptr, &pipeline_layout));

        auto vertShaderCode = readFile(expand_filename("vk/raycast.vert.spv"));
        auto fragShaderCode = readFile(expand_filename("vk/raycast.frag.spv"));

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
            sizeof(PackedAABB), // stride
            VK_VERTEX_INPUT_RATE_INSTANCE, // AABB array is instanced.
        };
        std::array<VkVertexInputAttributeDescription, 2> attributeDescriptions {
            VkVertexInputAttributeDescription {
              0, // location
              0, // binding
              VK_FORMAT_R32_SINT,
              offsetof(PackedAABB, packed_low) },
            VkVertexInputAttributeDescription {
              1, // location
              0, // binding
              VK_FORMAT_R32_SINT,
              offsetof(PackedAABB, packed_high) },
        };

        vertexInputInfo.vertexBindingDescriptionCount = 1;
        vertexInputInfo.vertexAttributeDescriptionCount = static_cast<uint32_t>(attributeDescriptions.size());
        vertexInputInfo.pVertexBindingDescriptions = &bindingDescription;
        vertexInputInfo.pVertexAttributeDescriptions = attributeDescriptions.data();

        VkPipelineInputAssemblyStateCreateInfo inputAssembly{};
        inputAssembly.sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO;
        inputAssembly.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP;
        inputAssembly.primitiveRestartEnable = VK_FALSE;

        VkViewport viewport{};
        viewport.x = 0.0f;
        viewport.y = 0.0f;
        viewport.width = 1.0f;
        viewport.height = 1.0f;
        viewport.minDepth = min_depth;
        viewport.maxDepth = max_depth;

        VkRect2D scissor{};
        scissor.offset = {0, 0};
        scissor.extent = VkExtent2D{ 1, 1 };

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
        pipelineInfo.flags = VK_DYNAMIC_STATE_VIEWPORT | VK_DYNAMIC_STATE_SCISSOR;
        pipelineInfo.stageCount = 2;
        pipelineInfo.pStages = shaderStages;
        pipelineInfo.pVertexInputState = &vertexInputInfo;
        pipelineInfo.pInputAssemblyState = &inputAssembly;
        pipelineInfo.pViewportState = &viewportState;
        pipelineInfo.pRasterizationState = &rasterizer;
        pipelineInfo.pMultisampleState = &multisampling;
        pipelineInfo.pDepthStencilState = &depthStencil;
        pipelineInfo.pColorBlendState = &colorBlending;
        pipelineInfo.layout = pipeline_layout;
        pipelineInfo.renderPass = render_pass;
        pipelineInfo.subpass = 0;
        pipelineInfo.basePipelineHandle = VK_NULL_HANDLE;

        NVVK_CHECK(vkCreateGraphicsPipelines(device, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &pipeline));

        vkDestroyShaderModule(device, fragShaderModule, nullptr);
        vkDestroyShaderModule(device, vertShaderModule, nullptr);
    }

    RaycastPipeline(RaycastPipeline&&) = delete;

    ~RaycastPipeline()
    {
        vkDestroyPipeline(device, pipeline, nullptr);
        vkDestroyPipelineLayout(device, pipeline_layout, nullptr);
        vkDestroyDescriptorSetLayout(device, set_layout, nullptr);
    }

    operator VkPipeline() const
    {
        return pipeline;
    }
};




// Dynamic viewport/scissor pipeline for drawing the background.
// No depth write or test, and the vertex shader hard-codes
// drawing a full-screen quad as a 4 vertex triangle strip.
struct BackgroundPipeline
{
    // We manage these.
    VkPipeline pipeline;
    VkPipelineLayout layout;

    // Borrowed pointer.
    VkDevice device;

    BackgroundPipeline(VkDevice dev_, VkRenderPass render_pass)
    {
        device = dev_;

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

        auto vertShaderCode = readFile(expand_filename("vk/background.vert.spv"));
        auto fragShaderCode = readFile(expand_filename("vk/background.frag.spv"));

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

        // Hard coded 4-vertex quad in vertex shader, so no vertex attributes.
        vertexInputInfo.vertexBindingDescriptionCount = 0;
        vertexInputInfo.vertexAttributeDescriptionCount = 0;
        vertexInputInfo.pVertexBindingDescriptions = nullptr;
        vertexInputInfo.pVertexAttributeDescriptions = nullptr;

        VkPipelineInputAssemblyStateCreateInfo inputAssembly{};
        inputAssembly.sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO;
        inputAssembly.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP;
        inputAssembly.primitiveRestartEnable = VK_FALSE;

        VkViewport viewport{};
        viewport.x = 0.0f;
        viewport.y = 0.0f;
        viewport.width = 1.0f;
        viewport.height = 1.0f;
        viewport.minDepth = min_depth;
        viewport.maxDepth = max_depth;

        VkRect2D scissor{};
        scissor.offset = {0, 0};
        scissor.extent = VkExtent2D{ 1, 1 };

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
        depthStencil.depthTestEnable = VK_FALSE;
        depthStencil.depthWriteEnable = VK_FALSE;
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
        pipelineInfo.flags = VK_DYNAMIC_STATE_VIEWPORT | VK_DYNAMIC_STATE_SCISSOR;
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

    BackgroundPipeline(BackgroundPipeline&&) = delete;

    ~BackgroundPipeline()
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
    // The order of these is very important (things need to be
    // destroyed in the right order).
    ScopedContext ctx;
    VkDevice dev;
    ScopedSurface glfw_surface;
    ScopedRenderPass render_pass;
    Framebuffers framebuffers;
    MeshPipeline mesh_pipeline;
    RaycastPipeline raycast_pipeline;
    BackgroundPipeline background_pipeline;

    // Last, as FrameManager destructor blocks the main queue,
    // allowing the aboves' destructors to run safely.
    FrameManager frame_manager;

    // Set in begin_frame.
    uint32_t width, height;

    // Graphics+Present queue used by the main thread, and its family.
    // TODO possible rename, as transfer_queue is now also used by
    // main thread only.
    VkQueue main_queue;
    uint32_t main_queue_family;

    // Transfer-only queue. Command buffers submitted to are recorded
    // and submitted on the main thread.
    VkQueue transfer_queue;
    uint32_t transfer_queue_family;

    // Command pool for allocating command buffers for transfer_queue.
    VkCommandPool transfer_cmd_pool;

    // Command buffer for submitting to transfer_queue, plus a fence
    // synchronizing access to it. This command buffer can be
    // re-recorded.
    struct TransferCmdBuffer
    {
        VkCommandBuffer buffer = VK_NULL_HANDLE;

        // Three possible states:
        //
        // non-nullptr, signalled: buffer submitted to queue,
        // retired. Command buffer may be re-recorded.
        //
        // non-nullptr, not-signalled: buffer submitted to queue, not
        // yet retired OR command buffer currently recording; this
        // fence will be signalled after buffer submitted and retired.
        //
        // nullptr: all other possibilities.
        std::shared_ptr<VkFence> p_fence;
    };

    // Pool of transfer-only command buffers. One is chosen for
    // current recording. Each frame the main thread will submit the
    // command buffer and choose a new command buffer from the array
    // for recording. NOTE: p_transfer_cmd_buffer may be nullptr if
    // all are in-flight (i.e. heavy load condition).
    std::array<TransferCmdBuffer, 24> transfer_cmd_buffer_array;
    TransferCmdBuffer* p_transfer_cmd_buffer = &transfer_cmd_buffer_array[0];

    // Number of chunk groups' RaycastStaging transfer commands
    // recorded to transfer command buffer, and max number recorded
    // before submitting.
    int32_t cmdbuf_group_transfer_count = 0;
    static constexpr int32_t max_cmdbuf_group_transfer_count = 7;

    // From FrameManager: primary command buffer that will be
    // submitted at frame-end, and the current swap chain image. Both
    // owned by the main thread.
    VkCommandBuffer frame_cmd_buffer = VK_NULL_HANDLE;
    nvvk::SwapChainImage current_swap_image;

    // From FrameManager: primary command buffer submitted before
    // frame_cmd_buffer.
    VkCommandBuffer pre_frame_buffer = VK_NULL_HANDLE;

    // Garbage produced by MeshStore/RaycastStore. Only destroyed
    // after the stores are destroyed, as only at that point will all
    // MeshEntry/RaycastEntry destructors be called (fully filling out
    // the list). This is ensured by RenderThread::destroy_stores.
    std::vector<MappedBuffer<MappedGroupMesh>>  unused_mesh_buffers;
    std::vector<MappedBuffer<MappedChunks>>     unused_chunk_buffers;
    std::vector<MappedBuffer<AABBs>>            unused_aabb_buffers;
    std::vector<std::pair<VkImage, VkImageView>>unused_voxel_images;
    std::vector<VkDescriptorPool>               descriptor_pools_to_free;
    std::vector<VkDescriptorSet>                unused_descriptor_sets;

    // VkDeviceMemory we're using that is to be freed in the destructor.
    std::vector<VkDeviceMemory> memory_to_free;

    // Push constant. Frame-constant fields set in begin_frame (too
    // small to be worth using a UBO, plus I'm lazy).
    PushConstant push_constant{};

    // Used in flush_transfer_cmd_buffers.
    bool cmd_buffer_fence_warning = false;
    bool p_transfer_cmd_buffer_nullptr_warning = false;
    size_t transfer_cmd_buffer_index = 0;

    RendererVk(RenderThread* thread, RenderArgs args) :
        RendererLogic<MeshEntry, MeshStaging, RaycastEntry, RaycastStaging>(
            thread,
            args),
        ctx(),
        dev(ctx.m_device),
        glfw_surface(
            ctx.m_instance,
            args.p_window->get_glfw_window()),
        render_pass(
            dev),
        framebuffers(
            ctx.m_physicalDevice,
            dev,
            render_pass),
        mesh_pipeline(
            dev,
            render_pass),
        raycast_pipeline(
            dev,
            render_pass),
        background_pipeline(
            dev,
            render_pass),
        frame_manager(
            ctx,
            glfw_surface.surface,
            glfw_surface.initial_width,
            glfw_surface.initial_height)
    {
        // As promised.
        thread_local_renderer = this;

        // Steal the queues from the ever-helpful nvvk::Context.
        // ctx.m_queueGCT is the same queue FrameManager is using.
        main_queue =        ctx.m_queueGCT;
        main_queue_family = ctx.m_queueGCT;
        transfer_queue =        ctx.m_queueT;
        transfer_queue_family = ctx.m_queueT;

        if (main_queue == VK_NULL_HANDLE) {
            throw std::runtime_error("Failed to find graphics+present queue");
        }
        if (transfer_queue == VK_NULL_HANDLE) {
            throw std::runtime_error("Failed to find transfer-only queue");
        }

        // Since in principle the pre_frame_buffer could have commands
        // recorded into it at any time (not just between begin and
        // end frame), we have to start recording right now.
        pre_frame_buffer = frame_manager.recordOneTimeCommandBuffer();

        // Initialize the command pool and command buffers for the
        // transfer-only queue.
        VkCommandPoolCreateInfo cmd_pool_info {
            VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO,
            nullptr,
            VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT,
            transfer_queue_family };
        NVVK_CHECK(vkCreateCommandPool(
            dev, &cmd_pool_info, nullptr, &transfer_cmd_pool));

        for (TransferCmdBuffer& wcb : transfer_cmd_buffer_array) {
            VkCommandBufferAllocateInfo info {
                VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,
                nullptr,
                transfer_cmd_pool,
                VK_COMMAND_BUFFER_LEVEL_PRIMARY,
                1 };
            NVVK_CHECK(vkAllocateCommandBuffers(dev, &info, &wcb.buffer));

            // wcb.p_fence == nullptr is correct.
            // Not recording or submitted to queue yet.
        }

        // Put one transfer command buffer in recording state.
        p_transfer_cmd_buffer = transfer_cmd_buffer_array.data();
        flush_transfer_cmd_buffer(true);
    }

    ~RendererVk()
    {
        // By the time this destructor runs, RendererBase has
        // destroyed MeshStore and RaycastStore, stopping all worker
        // threads.
        vkDestroyCommandPool(dev, transfer_cmd_pool, nullptr);

        for (auto mapped_buffer : unused_mesh_buffers) {
            vkDestroyBuffer(dev, mapped_buffer.buffer, nullptr);
        }
        for (auto mapped_buffer : unused_chunk_buffers) {
            vkDestroyBuffer(dev, mapped_buffer.buffer, nullptr);
        }
        for (auto mapped_buffer : unused_aabb_buffers) {
            vkDestroyBuffer(dev, mapped_buffer.buffer, nullptr);
        }

        for (std::pair<VkImage, VkImageView> pair : unused_voxel_images) {
            vkDestroyImageView(dev, pair.second, nullptr);
            vkDestroyImage(dev, pair.first, nullptr);
        }

        for (VkDescriptorPool pool : descriptor_pools_to_free) {
            vkDestroyDescriptorPool(dev, pool, nullptr);
        }

        for (VkDeviceMemory mem : memory_to_free) {
            vkFreeMemory(dev, mem, nullptr);
        }
    }

    // Stop recording and submit the transfer-only command buffer,
    // then pick a new command buffer, start recording it, and set it
    // as the command buffer to record transfer commands to.
    // first_time=true is used for initialization: skip the
    // end-record-and-submit in that case.
    //
    // Shoehorned hack: if we couldn't acquire a new command buffer
    // (all are in flight an their fences aren't yet signalled), don't
    // start recording a new command buffer, and set
    // p_transfer_cmd_buffer to nullptr. Subsequent calls to
    // flush_transfer_cmd_buffer can resolve this situation.
    void flush_transfer_cmd_buffer(bool first_time=false) noexcept
    {
        size_t& idx = transfer_cmd_buffer_index;

        // If we aren't stalled by a fence from a previous call, end
        // command buffer recording and select a new command
        // buffer. Since I was recording this buffer, the p_fence
        // should not be nullptr.
        if (!first_time and p_transfer_cmd_buffer != nullptr) {
            vkEndCommandBuffer(p_transfer_cmd_buffer->buffer);
            assert(p_transfer_cmd_buffer->p_fence != nullptr);

            // Submit to transfer-only queue.
            VkSubmitInfo submit {
                VK_STRUCTURE_TYPE_SUBMIT_INFO,
                nullptr,
                0, nullptr, nullptr,
                1, &p_transfer_cmd_buffer->buffer,
                0, nullptr };
            vkQueueSubmit(
                transfer_queue, 1, &submit, *p_transfer_cmd_buffer->p_fence);
            // fprintf(stderr, "Submitted transfer cmd buffer %i\n", idx);

            // Pick another command buffer for swap_in RaycastStaging
            // to record to.
            if (++idx >= transfer_cmd_buffer_array.size()) {
                idx = 0;
            }
        }
        p_transfer_cmd_buffer = &transfer_cmd_buffer_array[idx];

        // Wait on a fence, if neccessary. If we fail, null out
        // p_transfer_cmd_buffer, exit early, and hope the situation
        // fixes itself in a subsequent call.
        if (p_transfer_cmd_buffer->p_fence != nullptr) {
            VkFence* p_fence = p_transfer_cmd_buffer->p_fence.get();

            VkResult result = vkWaitForFences(dev, 1, p_fence, true, 0);
            if (result == VK_TIMEOUT) {
                if (!cmd_buffer_fence_warning) {
                    fprintf(stderr, "flush_transfer_cmd_buffer: no cmd buffer available\n");
                    cmd_buffer_fence_warning = true;
                }
                else {
                    assert(result == VK_SUCCESS);
                }
                p_transfer_cmd_buffer = nullptr;
                return;
            }
        }
        // else {
        //     fprintf(stderr, "Hello\n");
        // }

        // Now that recording will start, I need to make p_fence non-null.
        VkFenceCreateInfo fence_info {
            VK_STRUCTURE_TYPE_FENCE_CREATE_INFO,
            nullptr,
            0 };
        VkFence* p_fence = new VkFence;
        // static long fence_count; // Surprisingly useful
        // fence_count++;
        vkCreateFence(dev, &fence_info, nullptr, p_fence);
        auto del = [dev=dev] (VkFence* p_fence) {
            vkDestroyFence(dev, *p_fence, nullptr);
            // fprintf(stderr, "Fence count: %i\n", fence_count);
            // fence_count--;
        };
        p_transfer_cmd_buffer->p_fence.reset(p_fence, del);

        // Actually start recording said command buffer.
        VkCommandBufferBeginInfo begin_info {
            VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,
            nullptr,
            0,
            nullptr };
        vkBeginCommandBuffer(p_transfer_cmd_buffer->buffer, &begin_info);
    }

    void begin_frame() override
    {
        // Begin the frame, recreating the swap chain (hidden in
        // beginFrame) and framebuffers if needed. This converts the
        // GLFW window size to the actual framebuffer size (might
        // differ due to lag, etc.)
        int glfw_width, glfw_height;
        glfwGetFramebufferSize(glfw_surface.p_window, &glfw_width, &glfw_height);
        width = uint32_t(glfw_width), height = uint32_t(glfw_height);
        frame_manager.beginFrame(
            &frame_cmd_buffer,
            &current_swap_image,
            &width, &height);
        framebuffers.recreate_now_if_needed(frame_manager.getSwapChain());

        VkClearColorValue clear_color;
        clear_color.float32[0] = 0.0f;
        clear_color.float32[1] = 0.0f;
        clear_color.float32[2] = 0.0f;
        clear_color.float32[3] = 1.0f;
        std::array<VkClearValue, 2> clear_values{};
        clear_values[0].color = clear_color;            // Color attachment
        clear_values[1].depthStencil = { max_depth, 0 };// Depth attachment

        // Write the frame-constant push constant fields.
        push_constant.flags = 0;
        if (transforms.use_fog) {
            push_constant.flags |= MYRICUBE_FOG_BIT;
        }
        if (transforms.use_black_fog) {
            push_constant.flags |= MYRICUBE_BLACK_FOG_BIT;
        }
        if (transforms.chunk_debug) {
            push_constant.flags |= MYRICUBE_CHUNK_DEBUG_BIT;
        }

        push_constant.far_plane_squared =
            transforms.far_plane * transforms.far_plane;
        push_constant.raycast_thresh_squared =
            raycast_threshold * raycast_threshold;

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

        // Draw the background first.
        vkCmdBindPipeline(
            frame_cmd_buffer,
            VK_PIPELINE_BIND_POINT_GRAPHICS,
            background_pipeline);

        // Set viewport and scissors.
        VkViewport viewport {
            0, 0, float(width), float(height), min_depth, max_depth };
        vkCmdSetViewport(frame_cmd_buffer, 0, 1, &viewport);

        VkRect2D scissor { { 0, 0 }, { uint32_t(width), uint32_t(height) } };
        vkCmdSetScissor(frame_cmd_buffer, 0, 1, &scissor);

        // Set push constant matrix and eye vector, see
        // background.frag for real meanings.
        push_constant.mvp = glm::inverse(transforms.residue_vp_matrix);
        push_constant.eye_relative_group_origin =
            glm::vec4(transforms.eye_residue, 0.0f);
        vkCmdPushConstants(
            frame_cmd_buffer,
            background_pipeline.layout,
            VK_SHADER_STAGE_ALL_GRAPHICS,
            0, sizeof(PushConstant), &push_constant);

        // Draw full-screen quad.
        vkCmdDraw(frame_cmd_buffer, 4, 1, 0, 0);
    }


    /* MESH RENDERING IMPLEMENTATION */

    // Add to the primary command buffer commands to draw the list
    // of chunk groups with the instanced-voxel rendering method.
    // At this point, we are in the render pass begin_frame started.
    void draw_mesh_entries(
        const std::vector<std::pair<MeshEntry*, glm::ivec3>>& entries) override
    {
        glm::mat4 vp = transforms.residue_vp_matrix;
        glm::vec3 eye_residue = transforms.eye_residue;
        glm::ivec3 eye_group = transforms.eye_group;

        // Bind the mesh-drawing pipeline.
        vkCmdBindPipeline(
            frame_cmd_buffer,
            VK_PIPELINE_BIND_POINT_GRAPHICS,
            mesh_pipeline);

        // Set viewport and scissors.
        VkViewport viewport {
            0, 0, float(width), float(height), min_depth, max_depth };
        vkCmdSetViewport(frame_cmd_buffer, 0, 1, &viewport);

        VkRect2D scissor { { 0, 0 }, { uint32_t(width), uint32_t(height) } };
        vkCmdSetScissor(frame_cmd_buffer, 0, 1, &scissor);

        // Ready to draw the entries.
        for (auto pair : entries) {
            MeshEntry& entry = *pair.first;
            glm::ivec3 group_coord = pair.second;

            // Push the push constant. We need to update the MVP and
            // the eye position relative to this group.  The view
            // matrix only takes into account the eye's residue
            // coordinates, so the model position of the group
            // actually needs to be shifted by the eye's group coord.
            glm::vec3 model_offset = glm::vec3(group_coord - eye_group)
                                   * float(group_size);
            glm::mat4 m = glm::translate(glm::mat4(1.0f), model_offset);
            push_constant.mvp = vp * m;
            push_constant.eye_relative_group_origin = glm::vec4(
                eye_residue - model_offset, 0);

            vkCmdPushConstants(
                frame_cmd_buffer,
                mesh_pipeline.layout,
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

    // Implement the constructor and destructor for MeshEntry
    // (basically just a host-mapped buffer of MeshVoxelVertex).  My
    // plan is to allocate big blocks of buffers (block = shared
    // vulkan device memory) and only release the memory at the very
    // end. Hence, I need to recycle buffers (in unused_mesh_buffers)
    // to avoid an unbounded memory leak.

    // Construct the given MeshEntry. Only called by main thread.
    void constructor(MeshEntry* entry)
    {
      restart:
        // Recycle if possible.
        if (!unused_mesh_buffers.empty()) {
            entry->buffer = unused_mesh_buffers.back().buffer;
            entry->map = unused_mesh_buffers.back().map;
            unused_mesh_buffers.pop_back();
            return;
        }

        // Otherwise, we have to allocate a new block of buffers.
        constexpr size_t block_size = 80;
        memory_to_free.push_back(VK_NULL_HANDLE);
        unused_mesh_buffers.reserve(unused_mesh_buffers.size() + block_size);

        auto buffers = create_mapped_buffer_array<MappedGroupMesh>(
            ctx.m_physicalDevice,
            dev,
            sizeof(MappedGroupMesh),
            block_size,
            VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
            VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT
                | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
            &memory_to_free.back());

        for (auto buffer : buffers) {
            unused_mesh_buffers.push_back(buffer);
        }
        goto restart;
    }

    // Release resources of MeshEntry.  Only called by main thread.
    void destructor(MeshEntry* entry) noexcept
    {
        unused_mesh_buffers.push_back( { entry->buffer, entry->map } );
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
    // it can be swapped into the main cache. This is done by the main
    // thread's swap_in function in two steps, which are distinguished
    // by the command_recorded flag.
    //
    // First, record commands to the transfer-only queue to copy from
    // the staging buffer to the 3D VkImage
    //
    // Second, the main thread checks a VkFence to see if this is done
    // and swaps the VkImage into the main cache when ready.
    //
    // In both step, we also need to schedule a pipeline barrier in
    // order to do the image layout transition and queue ownership
    // transfer.

    // Draw the given list of chunk groups with the raycast-AABB
    // rendering method.
    void draw_raycast_entries(
        const std::vector<std::pair<RaycastEntry*, glm::ivec3>>& entries)
    override
    {
        glm::mat4 vp = transforms.residue_vp_matrix;
        glm::vec3 eye_residue = transforms.eye_residue;
        glm::ivec3 eye_group = transforms.eye_group;

        // Bind the raycast-drawing pipeline.
        vkCmdBindPipeline(
            frame_cmd_buffer,
            VK_PIPELINE_BIND_POINT_GRAPHICS,
            raycast_pipeline);

        // Set viewport and scissor.
        VkViewport viewport {
            0, 0, float(width), float(height), min_depth, max_depth };
        vkCmdSetViewport(frame_cmd_buffer, 0, 1, &viewport);

        VkRect2D scissor { { 0, 0 }, { uint32_t(width), uint32_t(height) } };
        vkCmdSetScissor(frame_cmd_buffer, 0, 1, &scissor);

        // Draw the raycast entries.
        for (auto pair : entries) {
            RaycastEntry& entry = *pair.first;
            glm::ivec3 group_coord = pair.second;

            // Push the push constant. We need to update the MVP and
            // eye position relative to this group.  The view matrix
            // only takes into account the eye's residue coordinates,
            // so the model position of the group actually needs to be
            // shifted by the eye's group coord.
            glm::vec3 model_offset = glm::vec3(group_coord - eye_group)
                                   * float(group_size);
            glm::mat4 m = glm::translate(glm::mat4(1.0f), model_offset);
            push_constant.mvp = vp * m;
            push_constant.eye_relative_group_origin = glm::vec4(
                eye_residue - model_offset, 0);

            vkCmdPushConstants(
                frame_cmd_buffer,
                raycast_pipeline.pipeline_layout,
                VK_SHADER_STAGE_ALL_GRAPHICS,
                0, sizeof(PushConstant), &push_constant);

            // Bind the instanced AABB buffer for this group in
            // binding point 0.
            VkDeviceSize offsets{0};
            vkCmdBindVertexBuffers(
                frame_cmd_buffer, 0, 1, &entry.aabb_buffer, &offsets);

            // Bind the chunk group voxels 3D image.
            vkCmdBindDescriptorSets(
                frame_cmd_buffer,
                VK_PIPELINE_BIND_POINT_GRAPHICS,
                raycast_pipeline.pipeline_layout,
                0, 1, &entry.voxels_descriptor,
                0, nullptr);

            // Draw the AABBs, there are as many as there are chunks
            // in a chunk group and each takes 14 vertices.
            uint32_t instances = edge_chunks * edge_chunks * edge_chunks;
            vkCmdDraw(frame_cmd_buffer, 14, instances, 0, 0);
        }
    }

    // Fill in the AABB buffer, 3D image, and descriptor set bound to
    // that image. We'll use the same block strategy as in MeshEntry.
    void constructor(RaycastEntry* entry)
    {
        // Init the AABB buffer.
      retry_buffer:
        // Recycle if possible.
        if (!unused_aabb_buffers.empty()) {
            entry->aabb_buffer = unused_aabb_buffers.back().buffer;
            entry->aabb_map = unused_aabb_buffers.back().map;
            unused_aabb_buffers.pop_back();
        }
        // Otherwise, we have to allocate a new block of buffers.
        else {
            constexpr size_t block_size = 4096;
            memory_to_free.push_back(VK_NULL_HANDLE);
            unused_aabb_buffers.reserve(
                unused_aabb_buffers.size() + block_size);
            auto buffers = create_mapped_buffer_array<AABBs>(
                ctx.m_physicalDevice,
                dev,
                sizeof(AABBs),
                block_size,
                VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
                VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT
                    | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                &memory_to_free.back());

            for (auto buffer : buffers) {
                unused_aabb_buffers.push_back(buffer);
            }
            goto retry_buffer;
        }

        // Init the 3D chunk group voxel image.
      retry_image:
        if (!unused_voxel_images.empty()) {
            entry->voxels_image = unused_voxel_images.back().first;
            entry->voxels_view = unused_voxel_images.back().second;
            unused_voxel_images.pop_back();
        }
        else {
            constexpr size_t block_size = 256;
            memory_to_free.push_back(VK_NULL_HANDLE);
            unused_voxel_images.reserve(
                unused_voxel_images.size() + block_size);
            auto images = create_chunk_group_images(
                ctx.m_physicalDevice,
                dev,
                block_size,
                chunk_group_voxels_image_format,
                VK_IMAGE_TILING_OPTIMAL,
                VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_STORAGE_BIT,
                VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                &memory_to_free.back());

            for (auto image : images) {
                unused_voxel_images.push_back(image);
            }
            goto retry_image;
        }

        // Init the descriptor set.
      retry_set:
        if (!unused_descriptor_sets.empty()) {
            entry->voxels_descriptor = unused_descriptor_sets.back();
            unused_descriptor_sets.pop_back();
        }
        else {
            auto pool_info = RaycastPipeline::descriptor_pool_info;
            auto block_size = RaycastPipeline::pool_sets;
            descriptor_pools_to_free.push_back(VK_NULL_HANDLE);
            VkDescriptorPool* p_pool = &descriptor_pools_to_free.back();
            unused_descriptor_sets.reserve(
                unused_descriptor_sets.size() + block_size);
            NVVK_CHECK(vkCreateDescriptorPool(
                dev, &pool_info, nullptr, p_pool));

            VkDescriptorSetAllocateInfo set_info {
                VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
                nullptr,
                *p_pool,
                1,
                &raycast_pipeline.set_layout };

            for (auto i = block_size; i != 0; --i) {
                VkDescriptorSet set;
                vkAllocateDescriptorSets(dev, &set_info, &set);
                unused_descriptor_sets.push_back(set);
            }
            goto retry_set;
        }

        // Bind descriptor set to point to the image.
        VkDescriptorImageInfo image_info {
            VK_NULL_HANDLE,
            entry->voxels_view,
            VK_IMAGE_LAYOUT_GENERAL };
        VkWriteDescriptorSet write {
            VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
            nullptr,
            entry->voxels_descriptor,
            0,
            0,
            1,
            VK_DESCRIPTOR_TYPE_STORAGE_IMAGE,
            &image_info,
            nullptr,
            nullptr };
        vkUpdateDescriptorSets(dev, 1, &write, 0, nullptr);
    }

    void destructor(RaycastEntry* entry) noexcept
    {
        unused_aabb_buffers.push_back(
            { entry->aabb_buffer, entry->aabb_map } );
        unused_voxel_images.push_back(
            { entry->voxels_image, entry->voxels_view } );
        unused_descriptor_sets.push_back(
            { entry->voxels_descriptor } );
    }

    void constructor(RaycastStaging* staging)
    {
        // Create the RaycastEntry (including sub-resources)
        staging->entry.reset(new RaycastEntry());

        // Create the memory-mapped staging buffer.
      retry_buffer:
        // Recycle if possible.
        if (!unused_chunk_buffers.empty()) {
            staging->chunk_buffer = unused_chunk_buffers.back().buffer;
            staging->chunk_map = unused_chunk_buffers.back().map;
            unused_chunk_buffers.pop_back();
            return;
        }

        // Otherwise, we have to allocate a new block of buffers.
        constexpr size_t block_size = 24;
        memory_to_free.push_back(VK_NULL_HANDLE);
        unused_chunk_buffers.reserve(unused_chunk_buffers.size() + block_size);

        auto buffers = create_mapped_buffer_array<MappedChunks>(
            ctx.m_physicalDevice,
            dev,
            sizeof(MappedChunks),
            block_size,
            VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
            VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT
                | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
            &memory_to_free.back());

        for (auto buffer : buffers) {
            unused_chunk_buffers.push_back(buffer);
        }
        goto retry_buffer;
    }

    void destructor(RaycastStaging* staging)
    {
        // Real destructor will take care of the RaycastEntry.

        // Recycle the staging buffer.
        unused_chunk_buffers.push_back(
            { staging->chunk_buffer, staging->chunk_map } );
    }

    // Fill in both the AABB and staging chunk buffers.
    void worker_stage(
        RaycastStaging* staging, const BinChunkGroup* group_ptr) override
    {
        // Fill in the AABB buffer and voxel staging buffer.
        auto& aabb_array = staging->entry->aabb_map->aabb_array;

        for (int zL = 0; zL < edge_chunks; ++zL) {
        for (int yL = 0; yL < edge_chunks; ++yL) {
        for (int xL = 0; xL < edge_chunks; ++xL) {
            // Load raw voxel data chunk-by-chunk onto the SSBO.
            const BinChunk& bin_chunk = group_ptr->chunk_array[zL][yL][xL];
            auto* p_source_chunk = &bin_chunk.voxel_array;
            auto* p_device_chunk =
                &staging->chunk_map->chunks[zL][yL][xL];
            auto sz = sizeof(BinChunk);
            assert(sz == sizeof(*p_source_chunk));
            assert(sz == sizeof(*p_device_chunk));
            memcpy(p_device_chunk, p_source_chunk, sz);

            // Also load the AABB for this chunk.
            auto chunk_residue = glm::ivec3(xL, yL, zL) * chunk_size;
            PackedAABB aabb(bin_chunk, chunk_residue);
            aabb_array[zL][yL][xL] = aabb;
        }
        }
        }
    }

    // If not yet done, record a command to the transfer queue to copy
    // from the staging buffer to the 3D voxel image. Otherwise, check
    // to see if this transfer has finished, and swap into the main
    // cache if so (doing a queue ownership transfer from the transfer
    // queue to the main queue).
    bool swap_in(RaycastStaging* staging,
                 std::unique_ptr<RaycastEntry>* p_uptr_entry) override
    {
        if (staging->command_recorded) {
            return swap_in_if_ready(staging, p_uptr_entry);
        }
        else {
            record_cmd(staging);
            return false;
        }
    }

    // Possibility 1 of swap_in RaycastStaging.
    void record_cmd(RaycastStaging* staging)
    {
        if (p_transfer_cmd_buffer == nullptr) {
            if (!p_transfer_cmd_buffer_nullptr_warning) {
                fprintf(stderr, "record_cmd(RaycastStaging*): slowed by nullptr\n");
                p_transfer_cmd_buffer_nullptr_warning = true;
            }
            return;
        }

        RaycastEntry& staged_entry = *staging->entry;

        // Transition chunk group voxels 3D image to transfer dst layout.
        VkImageMemoryBarrier barrier {
            VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER,
            nullptr,
            VK_ACCESS_MEMORY_WRITE_BIT,
            VK_ACCESS_TRANSFER_WRITE_BIT,
            VK_IMAGE_LAYOUT_UNDEFINED,
            VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
            VK_QUEUE_FAMILY_IGNORED,
            VK_QUEUE_FAMILY_IGNORED,
            staged_entry.voxels_image,
            color_range };
        vkCmdPipelineBarrier(
            p_transfer_cmd_buffer->buffer,
            VK_PIPELINE_STAGE_ALL_COMMANDS_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            0,
            0, nullptr,
            0, nullptr,
            1, &barrier);

        // Copy chunk-by-chunk from staging buffer to image.
        for (int z = 0; z < edge_chunks; ++z) {
        for (int y = 0; y < edge_chunks; ++y) {
        for (int x = 0; x < edge_chunks; ++x) {
            VkBufferImageCopy region;
            region.bufferOffset =
                reinterpret_cast<char*>(&staging->chunk_map->chunks[z][y][x])
              - reinterpret_cast<char*>(&staging->chunk_map->chunks[0][0][0]);
            region.bufferRowLength = chunk_size;
            region.bufferImageHeight = chunk_size;
            region.imageSubresource = { VK_IMAGE_ASPECT_COLOR_BIT, 0, 0, 1 };
            region.imageOffset = { x*chunk_size, y*chunk_size, z*chunk_size };
            region.imageExtent = { chunk_size, chunk_size, chunk_size };

            vkCmdCopyBufferToImage(
                p_transfer_cmd_buffer->buffer,
                staging->chunk_buffer,
                staged_entry.voxels_image,
                VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                1,
                &region);
        }
        }
        }

        // Transition layout to general layout and transfer ownership
        // to the main queue family (used for graphics and swap chain
        // present).
        barrier = VkImageMemoryBarrier {
            VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER,
            nullptr,
            VK_ACCESS_TRANSFER_WRITE_BIT,
            VK_ACCESS_MEMORY_READ_BIT,
            VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
            VK_IMAGE_LAYOUT_GENERAL,
            transfer_queue_family,
            main_queue_family,
            staged_entry.voxels_image,
            color_range };
        vkCmdPipelineBarrier(
            p_transfer_cmd_buffer->buffer,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT,
            0,
            0, nullptr,
            0, nullptr,
            1, &barrier);

        // *p_transfer_cmd_buffer->p_fence is signalled when the command
        // buffer we are recording to is submitted and retired, so
        // make swapping into the main cache conditional on said fence
        // being signalled.
        staging->p_fence = p_transfer_cmd_buffer->p_fence;

        // Indicate cmd has been recorded.
        staging->command_recorded = true;
        // fprintf(stderr, "Recorded transfer command %i %p\n", cmdbuf_group_transfer_count, staging);

        // Submit the command buffer if it's full.
        if (++cmdbuf_group_transfer_count >= max_cmdbuf_group_transfer_count) {
            flush_transfer_cmd_buffer();
            cmdbuf_group_transfer_count = 0;
        }
    }

    // Possibility 2 of swap_in RaycastStaging.
    bool swap_in_if_ready(RaycastStaging* staging,
                          std::unique_ptr<RaycastEntry>* p_uptr_entry)
    {
        // Check if the entry is actually ready to be swapped in.
        auto result = vkWaitForFences(dev, 1, staging->p_fence.get(), true, 0);
        if (VK_TIMEOUT == result) {
            return false;
        }
        else {
            assert(VK_SUCCESS == result);
        }

        // Reset state of the RaycastStaging.
        staging->p_fence.reset();
        staging->command_recorded = false;

        RaycastEntry& staged_entry = *staging->entry;

        // Transfer ownership from the transfer queue to the main queue.
        VkImageMemoryBarrier barrier = VkImageMemoryBarrier {
            VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER,
            nullptr,
            VK_ACCESS_TRANSFER_WRITE_BIT,
            VK_ACCESS_MEMORY_READ_BIT,
            VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
            VK_IMAGE_LAYOUT_GENERAL,
            transfer_queue_family,
            main_queue_family,
            staged_entry.voxels_image,
            color_range };
        vkCmdPipelineBarrier(
            pre_frame_buffer,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VK_PIPELINE_STAGE_ALL_GRAPHICS_BIT,
            0,
            0, nullptr,
            0, nullptr,
            1, &barrier);

        // Swap into the main cache.
        std::swap(staging->entry, *p_uptr_entry);
        if (staging->entry == nullptr) {
            staging->entry.reset(new RaycastEntry());
        }
        return true;
    }



    void end_frame() override
    {
        // Submit any unsubmitted RaycastStaging transfer commands.
        // TODO I have to do this anyway even if 0 commands were
        // recorded, as I also use the flush function to check the
        // fence status (this is important if all available cmd
        // buffers are in flight). I may want to fix this
        // inefficiency.
        flush_transfer_cmd_buffer();

        // End the before-frame command buffer and submit, then schedule
        // it to be freed one frame later.
        vkEndCommandBuffer(pre_frame_buffer);
        VkSubmitInfo submit = {
            VK_STRUCTURE_TYPE_SUBMIT_INFO,
            nullptr,
            0, nullptr, nullptr,
            1, &pre_frame_buffer,
            0, 0 };
        NVVK_CHECK(vkQueueSubmit(main_queue, 1, &submit, VK_NULL_HANDLE));

        frame_manager.addFrameGarbage(
        [buffer=pre_frame_buffer] (const FrameManager& m) {
            vkFreeCommandBuffers(m.getDevice(), m.getCommandPool(), 1, &buffer);
        } );

        // End the one render pass there is.
        vkCmdEndRenderPass(frame_cmd_buffer);

        // Ends cmd buffer recording, submits it, and presents the
        // rendered image.
        frame_manager.endFrame(frame_cmd_buffer);
        frame_cmd_buffer = VK_NULL_HANDLE;

        // Get a new pre-frame command buffer.
        pre_frame_buffer = frame_manager.recordOneTimeCommandBuffer();
    }

    void main_thread_wait_idle() override
    {
        vkQueueWaitIdle(main_queue);
    }

    void wait_idle() override
    {
        vkQueueWaitIdle(main_queue);
        vkQueueWaitIdle(transfer_queue);
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

void constructor(MeshEntry* entry)
{
    thread_local_renderer->constructor(entry);
}

void destructor(MeshEntry* entry)
{
    thread_local_renderer->destructor(entry);
}

void constructor(RaycastEntry* entry)
{
    thread_local_renderer->constructor(entry);
}

void destructor(RaycastEntry* entry)
{
    thread_local_renderer->destructor(entry);
}

void constructor(RaycastStaging* staging)
{
    thread_local_renderer->constructor(staging);
}

void destructor(RaycastStaging* staging)
{
    thread_local_renderer->destructor(staging);
}

} // end anonymous namespace
