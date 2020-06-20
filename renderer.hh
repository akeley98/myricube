// Functions for drawing the voxel world. As typical for OpenGL, these
// are NOT threadsafe, even for drawing different worlds in parallel.

#ifndef MYRICUBE_RENDERER_HH_
#define MYRICUBE_RENDERER_HH_

namespace myricube {

class RaycastStore;
class MeshStore;

RaycastStore* new_raycast_store();
MeshStore* new_mesh_store();
void delete_raycast_store(RaycastStore*);
void delete_mesh_store(MeshStore*);

class Camera;
class VoxelWorld;

// Render, to the current framebuffer, chunks near the camera using
// the conventional mesh-based algorithm.
void render_world_mesh_step(const VoxelWorld&, Camera&);

// Render, to the current framebuffer, chunks around the camera
// using the AABB-raycast algorithm. Chunks that are near
// enough to have been drawn using the mesh algorithm will
// not be re-drawn.
void render_world_raycast_step(const VoxelWorld&, Camera&);

} // end namespace
#endif /* !MYRICUBE_RENDERER_HH_ */
