#include "myricube.hh"

#include <stdio.h>
#include <utility>
#include <vector>

#include "chunk.hh"
#include "renderer.hh"
#include "glad/glad.h"

namespace myricube {

static constexpr int edge_chunks = group_size / chunk_size;

// Single vertex of the VBO for a chunk's mesh.
struct MeshVertex
{
    uint32_t packed_vertex;
    uint32_t packed_color;

    // Pack this vertex given the vertex's residue coordinates
    // (i.e. local coordinate in coordinate system of a chunk group)
    // and 24-bit color.
    //
    // This constructor defines the packed layout.
    MeshVertex(uint8_t x, uint8_t y, uint8_t z,
               uint8_t red, uint8_t green, uint8_t blue)
    {
        static_assert(group_size <= 255,
                      "group too big for 8-bit unsigned coordinates.");
        packed_color = x | uint32_t(y) << 8 | uint32_t(z) << 16;
        packed_vertex = blue | uint32_t(green) << 8 | uint32_t(red) << 16;
    }
};

// CPU-side copy of the mesh for a single chunk.
struct ChunkMesh
{
    // CPU vertex list.
    std::vector<MeshVertex> verts;

    // Offset of first byte in the VBO where the GPU copy of the
    // vertex list is stored.
    GLsizeiptr vbo_byte_offset;

    // Number of bytes reserved in the VBO for this chunk.
    GLsizeiptr vbo_bytes_reserved;
};

// Handle for GPU resources for the mesh of one chunk group.  The
// CPU-side state of this MeshEntry should always accurately describe
// the state of the copy of this data on the GPU and GL. (i.e. if you
// modify the data within, it is your responsibility to update the GPU
// state at the same time).
struct MeshEntry
{
    // World ID and group coordinate of the chunk group this mesh is
    // for. This can be used to tell if this entry is correct or needs
    // to be replaced in MeshStore.
    uint64_t world_id = uint64_t(-1);
    glm::ivec3 group_coord;

    // mesh_array[z][y][x] is the mesh for ChunkGroup::chunk_array[z][y][x].
    ChunkMesh mesh_array[edge_chunks][edge_chunks][edge_chunks];

    // OpenGL name for the VBO used to store the meshes of this chunk group.
    // 0 when not yet allocated.
    GLuint vbo_name = 0;

    // Size in bytes of the vbo's data store on the GPU.
    GLsizeiptr vbo_bytes = 0;

    // Bytes used for the VBO (i.e. the byte offset where new data can be
    // written to the vbo).
    GLsizeiptr bytes_used = 0;

    MeshEntry() = default;

    // Get some RAII going.
    ~MeshEntry()
    {
        glDeleteBuffers(1, &vbo_name);
        vbo_name = 0;
        vbo_bytes = 0;
        bytes_used = 0;
    }

    // Disable moves, but, see swap below. (Harkens back to C++98 hacks :P )
    MeshEntry(MeshEntry&&) = delete;
};

void swap(MeshEntry& left, MeshEntry& right)
{
    using std::swap;
    swap(left.world_id, right.world_id);
    swap(left.group_coord, right.group_coord);
    swap(left.mesh_array, right.mesh_array);
    swap(left.vbo_name, right.vbo_name);
    swap(left.vbo_bytes, right.vbo_bytes);
    swap(left.bytes_used, right.bytes_used);
}

// Cache of MeshEntry/RaycastEntry objects. The idea is to store entry
// structs for chunk groups near the camera. They are stored in a
// wraparound fashion (i.e.  using modular arithmetic on group
// coordinates) so that as the camera moves, the new groups being
// rendered overwrites the old groups no longer being rendered.
template <class Entry, int N>
class BaseStore
{
    // N should be made tunable later depending on
    // Camera::raycast_threshold or Camera::far_plane.

    Entry entry_array[N][N][N];

  public:
    // Return a pointer to the location that the Entry for the chunk
    // group with the given group coordinate would be, IF it exists.
    // (i.e. that location may be blank or occupied by an entry for
    // another group).
    Entry* cached_location(glm::ivec3 group_coord)
    {
        uint32_t x = (uint32_t(group_coord.x) ^ 0x1000'0000) % N;
        uint32_t y = (uint32_t(group_coord.y) ^ 0x1000'0000) % N;
        uint32_t z = (uint32_t(group_coord.z) ^ 0x1000'0000) % N;
        return &entry_array[z][y][x];
    }

    // Return a valid, updated Entry for the given chunk group
    // extracted from the given world. TODO document templates needed
    // for this to work.
    template <typename EntryFiller>
    Entry* update(PositionedChunkGroup& pcg, VoxelWorld& world)
    {
        glm::ivec3 gc = group_coord(pcg);
        Entry* entry = cached_location(gc);

        // If the Entry in the array already corresponds to the given
        // chunk group; just update it (deal with dirty chunks) and
        // return.
        if (entry->world_id == world.id()
        and entry->group_coord == group_coord(pcg)) {
            EntryFiller::update(pcg, world, entry);
            return entry;
        }

        // Otherwise, we need to evict the current entry and repopulate it.
        fprintf(stderr, "Evicting MeshStore entry for chunk group "
                        "(%i %i %i) of world %li\n",
                        int(entry->group_coordinate.x),
                        int(entry->group_coordinate.y),
                        int(entry->group_coordinate.z),
                        long(entry->world_id));

        EntryFiller::replace(pcg, world, entry);
        return entry;
    }
};

class MeshStore : BaseStore<MeshEntry, 16> { };

MeshStore* new_mesh_store()
{
    return new MeshStore;
}

void delete_mesh_store(MeshStore* mesh_store)
{
    delete mesh_store;
}

// Renderer class. Instantiate it with the camera and world to render
// and use it once.
class Renderer
{
    Camera& camera;
    VoxelWorld& world;
  public:
    Renderer(Camera& camera_, VoxelWorld& world_) :
        camera(camera_), world(world_) { }
  private:
    // Given a chunk, make the mesh needed to render it. I also need
    // the residue coordinate of the lower-left of the chunk (i.e.
    // the position of the chunk relative to the chunk group it is
    // in) because the mesh is defined to be in residue coordinates.
    //
    // This is just ugly.
    //
    // TODO: consider properly handling edge cases with other chunks?
    // It's less easy than it seems, as the mesh can change not only
    // due to changes in this chunk, but also changes in neighbor
    // chunks that reveal previously hidden voxel faces. This may
    // not be worth it, so for now, I don't include a VoxelWorld
    // argument needed to make this happen.
    static std::vector<MeshVertex> make_mesh_verts(const Chunk& chunk,
                                                   glm::ivec3 chunk_residue)
    {
        // Look up whether the voxel at the given coordinate
        // (relative to the lower-left of this chunk) is visible.
        // Act as if voxels outside the chunk are always invisible.
        auto visible_block = [&] (glm::ivec3 coord) -> bool
        {
            if (coord.x < 0 or coord.x >= chunk_size
             or coord.y < 0 or coord.y >= chunk_size
             or coord.z < 0 or coord.z >= chunk_size) return false;
            // Note: the masking in Chunk::operator () won't mess up
            // this coord. (I'm kind of violating my own comment in
            // Chunk::operator() because coord won't actually be in
            // the chunk unless that chunk is at (0,0,0)).
            return chunk(coord).visible;
        };

        MeshVertex reference(
                chunk_residue.x,
                chunk_residue.y,
                chunk_residue.z,
                0, 0, 0);
        std::vector<MeshVertex> verts;

        auto verts_add_block = [&] (glm::ivec3 coord)
        {
            Voxel v = chunk(coord);
            if (!v.visible) return;

            auto r = v.red;
            auto g = v.green;
            auto b = v.blue;
            uint8_t x = uint8_t(coord.x) + chunk_residue.x;
            uint8_t y = uint8_t(coord.y) + chunk_residue.y;
            uint8_t z = uint8_t(coord.z) + chunk_residue.z;
            uint8_t x1 = x+1;
            uint8_t y1 = y+1;
            uint8_t z1 = z+1;

            if (!visible_block(coord + glm::ivec3(-1, 0, 0))) {
                verts.emplace_back(x,  y1, z1, r, g, b);
                verts.emplace_back(x,  y1, z,  r, g, b);
                verts.emplace_back(x,  y,  z1, r, g, b);
                verts.emplace_back(x,  y1, z,  r, g, b);
                verts.emplace_back(x,  y,  z,  r, g, b);
                verts.emplace_back(x,  y,  z1, r, g, b);
            }
            if (!visible_block(coord + glm::ivec3(1, 0, 0))) {
                verts.emplace_back(x1, y,  z,  r, g, b);
                verts.emplace_back(x1, y1, z,  r, g, b);
                verts.emplace_back(x1, y,  z1, r, g, b);
                verts.emplace_back(x1, y1, z,  r, g, b);
                verts.emplace_back(x1, y1, z1, r, g, b);
                verts.emplace_back(x1, y,  z1, r, g, b);
            }
            if (!visible_block(coord + glm::ivec3(0, -1, 0))) {
                verts.emplace_back(x,  y,  z,  r, g, b);
                verts.emplace_back(x1, y,  z1, r, g, b);
                verts.emplace_back(x,  y,  z1, r, g, b);
                verts.emplace_back(x,  y,  z,  r, g, b);
                verts.emplace_back(x1, y,  z,  r, g, b);
                verts.emplace_back(x1, y,  z1, r, g, b);
            }
            if (!visible_block(coord + glm::ivec3(0, 1, 0))) {
                verts.emplace_back(x1, y1, z1, r, g, b);
                verts.emplace_back(x1, y1, z,  r, g, b);
                verts.emplace_back(x,  y1, z,  r, g, b);
                verts.emplace_back(x,  y1, z1, r, g, b);
                verts.emplace_back(x1, y1, z1, r, g, b);
                verts.emplace_back(x,  y1, z,  r, g, b);
            }
            if (!visible_block(coord + glm::ivec3(0, 0, -1))) {
                verts.emplace_back(x,  y,  z,  r, g, b);
                verts.emplace_back(x,  y1, z,  r, g, b);
                verts.emplace_back(x1, y,  z,  r, g, b);
                verts.emplace_back(x1, y1, z,  r, g, b);
                verts.emplace_back(x1, y,  z,  r, g, b);
                verts.emplace_back(x,  y1, z,  r, g, b);
            }
            if (!visible_block(coord + glm::ivec3(0, 0, 1))) {
                verts.emplace_back(x,  y1, z1, r, g, b);
                verts.emplace_back(x1, y,  z1, r, g, b);
                verts.emplace_back(x1, y1, z1, r, g, b);
                verts.emplace_back(x1, y,  z1, r, g, b);
                verts.emplace_back(x,  y1, z1, r, g, b);
                verts.emplace_back(x,  y,  z1, r, g, b);
            }
        };

        for (int z = 0; z < chunk_size; ++z) {
            for (int y = 0; y < chunk_size; ++y) {
                for (int x = 0; x < chunk_size; ++x) {
                    verts_add_block(glm::ivec3(x, y, z));
                }
            }
        }
        return verts;
    }

    // Return the number of bytes that will be reserved in a VBO to
    // store the given array of mesh vertices. This is not the number
    // needed, but the number that will be allocated when we have a
    // choice.
    static GLsizeiptr verts_reserved_bytes(
        const std::vector<MeshVertex>& verts)
    {
        return GLsizeiptr((128+verts.size()) * sizeof(verts[0]));
    }

    // Return the actual number of bytes needed to store the given
    // mesh vertex array.
    static GLsizeiptr verts_needed_bytes(
        const std::vector<MeshVertex>& verts)
    {
        return GLsizeiptr(verts.size() * sizeof(verts[0]));
    }

    // Fill the given MeshEntry with mesh data for the given chunk
    // group from the given world.
    static void replace(PositionedChunkGroup& pcg,
                        VoxelWorld& world,
                        MeshEntry* entry)
    {
        entry->world_id = world.id();
        entry->group_coord = group_coord(pcg);

        GLsizeiptr bytes_to_alloc = 0;
        for (int z = 0; z < edge_chunks; ++z) {
            for (int y = 0; y < edge_chunks; ++y) {
                for (int x = 0; x < edge_chunks; ++x) {
                    ChunkMesh& mesh = entry->mesh_array[z][y][x];
                    Chunk& chunk = group(pcg).chunk_array[z][y][x];
                    glm::ivec3 residue = glm::ivec3(x,y,z) * chunk_size;
                    std::vector<MeshVertex> verts =
                        make_mesh_verts(chunk, residue);
                    chunk.mesh_dirty = false;
                    auto reserved_bytes = verts_reserved_bytes(verts);

                    mesh.vbo_byte_offset = bytes_to_alloc;
                    bytes_to_alloc += reserved_bytes;
                    mesh.vbo_bytes_reserved = reserved_bytes;
                    mesh.verts = std::move(verts);
                }
            }
        }

        if (entry->vbo_name == 0 or bytes_to_alloc > entry->vbo_bytes) {
            if (entry->vbo_name == 0) glGenBuffers(1, &entry->vbo_name);
            glBindBuffer(GL_ARRAY_BUFFER, entry->vbo_name);
            glBufferData(GL_ARRAY_BUFFER, bytes_to_alloc,
                         nullptr, GL_DYNAMIC_DRAW);
            entry->vbo_bytes = bytes_to_alloc;
            entry->bytes_used = bytes_to_alloc;
            PANIC_IF_GL_ERROR;
        }

        // Parasitically depend on the MeshEntry update function.
        // This is a bit dicey -- the main violation is that MeshEntry
        // doesn't accurately describe GPU state now, as the actual
        // data has not been uploaded yet.
        update(pcg, world, entry, true);
    }

    // Update the given mesh entry with new data from dirty chunks.
    // (The always_dirty flag is used to allow replace to depend on
    // this function, if true, re-upload all chunks' data).
    //
    // Requires that the MeshEntry corresponds to the given chunk
    // group of the given world.
    static void update(PositionedChunkGroup& pcg,
                       VoxelWorld& world,
                       MeshEntry* entry,
                       bool always_dirty = false) noexcept
    {
        assert(entry->world_id == world.id());
        assert(entry->group_coord == group_coord(pcg));
        assert(entry->vbo_name != 0);

        // First step is to recompute the mesh vertices of any chunk
        // that needs it. This is also when we find out if a
        // reallocation needs to be done due to some mesh growing too
        // much.
        //
        // Don't reset the chunk dirty flag yet since we'll need it to
        // know if re-upload has to happen.
        bool requires_reallocation = false;
        GLsizeiptr bytes_to_alloc = 0; // Only used in case of realloc.
        for (int z = 0; z < edge_chunks; ++z) {
            for (int y = 0; y < edge_chunks; ++y) {
                for (int x = 0; x < edge_chunks; ++x) {
                    ChunkMesh& mesh = entry->mesh_array[z][y][x];
                    std::vector<MeshVertex>& verts = mesh.verts;
                    Chunk& chunk = group(pcg).chunk_array[z][y][x];
                    glm::ivec3 residue = glm::ivec3(x,y,z) * chunk_size;
                    bool dirty = always_dirty | chunk.mesh_dirty;

                    if (dirty)  {
                        verts = make_mesh_verts(chunk, residue);
                    }
                    auto verts_bytes = verts_needed_bytes(verts);
                    requires_reallocation |=
                        (verts_bytes > mesh.vbo_bytes_reserved);
                    bytes_to_alloc += verts_reserved_bytes(verts);
                }
            }
        }

        // Now actually upload the data. This depends on whether
        // reallocation is needed or not. We also clear the mesh dirty
        // flags at this time, now that we're actually uploading the
        // mesh data.
        if (requires_reallocation) {
            // If realloc is required, resize the VBO and re-upload
            // every chunk's mesh (I may find a more efficient way one
            // day).
            glBindBuffer(GL_ARRAY_BUFFER, entry->vbo_name);
            glBufferData(GL_ARRAY_BUFFER, bytes_to_alloc,
                         nullptr, GL_DYNAMIC_DRAW);
            entry->vbo_bytes = bytes_to_alloc;
            entry->bytes_used = bytes_to_alloc;
            PANIC_IF_GL_ERROR;

            GLsizeiptr offset = 0;

            for (int z = 0; z < edge_chunks; ++z) {
                for (int y = 0; y < edge_chunks; ++y) {
                    for (int x = 0; x < edge_chunks; ++x) {
                        ChunkMesh& mesh = entry->mesh_array[z][y][x];
                        std::vector<MeshVertex>& verts = mesh.verts;
                        Chunk& chunk = group(pcg).chunk_array[z][y][x];
                        chunk.mesh_dirty = false;

                        auto bytes_reserved = verts_reserved_bytes(verts);
                        auto bytes_needed = verts_needed_bytes(verts);
                        mesh.vbo_byte_offset = offset;
                        mesh.vbo_bytes_reserved = bytes_reserved;

                        glBufferSubData(GL_ARRAY_BUFFER, offset,
                                        bytes_needed, verts.data());
                        offset += bytes_reserved;
                    }
                }
            }
            PANIC_IF_GL_ERROR;
        }
        // No reallocation case. Just re-upload the dirty chunks,
        // except when always_dirty is true.
        else {
            glBindBuffer(GL_ARRAY_BUFFER, entry->vbo_name);
            PANIC_IF_GL_ERROR;
            for (int z = 0; z < edge_chunks; ++z) {
                for (int y = 0; y < edge_chunks; ++y) {
                    for (int x = 0; x < edge_chunks; ++x) {
                        ChunkMesh& mesh = entry->mesh_array[z][y][x];
                        std::vector<MeshVertex>& verts = mesh.verts;
                        Chunk& chunk = group(pcg).chunk_array[z][y][x];
                        bool dirty = always_dirty | chunk.mesh_dirty;

                        if (!dirty) continue;

                        chunk.mesh_dirty = false;
                        auto bytes_needed = verts_needed_bytes(verts);
                        assert(mesh.vbo_byte_offset + bytes_needed <= entry->vbo_bytes);
                        assert(bytes_needed <= mesh.vbo_bytes_reserved);
                        glBufferSubData(GL_ARRAY_BUFFER,
                                        mesh.vbo_byte_offset,
                                        bytes_needed,
                                        verts.data());
                    }
                }
            }
            PANIC_IF_GL_ERROR;
        }
    }
};

} // end namespace
