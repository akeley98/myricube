// Camera class. Mostly just some vectors and matrices and stuff;
// however, I also associate the GPU resource management classes with
// the camera. Since the set of GPU-loaded chunks depends on camera
// position, this seems like a somewhat natural place to store it.

#ifndef MYRICUBE_CAMERA_HH_
#define MYRICUBE_CAMERA_HH_

#include "myricube.hh"

#include <assert.h>
#include <math.h>
#include <stddef.h>

#include "renderer.hh"

namespace myricube {

class Camera
{
    // Mysterious GPU storage managers along for the ride.
    RaycastStore* ptr_raycast_store = nullptr;
    MeshStore* ptr_mesh_store = nullptr;

    // Position of the eye.
    glm::dvec3 eye = glm::dvec3(0.0, 0.0, 0.0);

    // Near and far plane for z-depth.
    // Far plane influences the chunk render distance.
    float near_plane = 0.1f;
    int far_plane = 384;

    // (roughly) minimum distance from the camera that a chunk needs
    // to be to switch from mesh to raycast graphics.
    // Keep as int to avoid rounding errors in distance culling.
    int raycast_threshold = 120;

    // TODO: I've added setters for the near/raycast/far planes, but
    // to do this safely I need to make MeshStore::N and
    // RaycastStore::N configurable.

    // Horizontal and vertical angle camera is pointed in.
    float theta = 1.5707f, phi = 1.5707f;

    // Field of view (y direction), radians.
    float fovy_radians = 1.0f;

    // Window size in pixels.
    int window_x = 1, window_y = 1;

    // Fog setting.
    bool fog_enabled = true;

    // *** True when members below need to be recomputed due to ***
    // *** changes in members above.                            ***
    bool dirty = true;

    // Group and residue coordinates of the eye.
    glm::ivec3 eye_group;
    glm::vec3 eye_residue;

    // Frenet frame of the camera (unit vectors).
    glm::vec3 forward_normal_vector;
    glm::vec3 right_vector;
    glm::vec3 up_vector;

    // View matrix: this ONLY accounts for eye_residue vector, not
    // the eye_group vector. i.e. every chunk needs to be manually
    // shifted to account for the eye_group when rendered. This is
    // to avoid catastrophic cancellation.
    glm::mat4 residue_view_matrix;

    // Projection matrix.
    glm::mat4 projection_matrix;

    // Projection * View matrix
    glm::mat4 residue_vp_matrix;

  public:
    // Maximum number of 3D voxel textures that may be sent to the GPU
    // per frame. Mystery parameter.
    size_t max_raycast_evict = 5;

    // Respond if needed to dirty flag and recompute derived data.
    void fix_dirty()
    {
        if (!dirty) return;
        split_coordinate(eye, &eye_group, &eye_residue);

        forward_normal_vector = glm::vec3(
            sinf(phi) * cosf(theta),
            cosf(phi),
            sinf(phi) * sinf(theta));

        right_vector = glm::normalize(
            glm::cross(forward_normal_vector, glm::vec3(0, 1, 0)));
        up_vector = glm::cross(right_vector, forward_normal_vector);

        residue_view_matrix = glm::lookAt(eye_residue,
                                          eye_residue + forward_normal_vector,
                                          glm::vec3(0, 1, 0));

        projection_matrix = glm::perspective(
            float(fovy_radians),
            float(window_x) / window_y,
            float(near_plane),
            float(far_plane));

        residue_vp_matrix = projection_matrix * residue_view_matrix;

        dirty = false;
    }

    Camera() = default;

    ~Camera()
    {
        unload_gpu_storage();
    }

    Camera(Camera&&) = delete;

    // GPU Storage functions.
    RaycastStore& get_raycast_store()
    {
        if (ptr_raycast_store == nullptr) {
            ptr_raycast_store = new_raycast_store();
        }
        return *ptr_raycast_store;
    }

    MeshStore& get_mesh_store()
    {
        if (ptr_mesh_store == nullptr) {
            ptr_mesh_store = new_mesh_store();
        }
        return *ptr_mesh_store;
    }

    void unload_gpu_storage()
    {
        delete_raycast_store(ptr_raycast_store);
        ptr_raycast_store = nullptr;
        delete_mesh_store(ptr_mesh_store);
        ptr_mesh_store = nullptr;
    }

    // Getters and setters for user-specified camera data.
    glm::dvec3 get_eye() const
    {
        return eye;
    }

    float get_near_plane() const
    {
        return near_plane;
    }

    void set_near_plane(float in)
    {
        near_plane = in;
        dirty = true;
    }

    int get_far_plane() const
    {
        return far_plane;
    }

    void set_far_plane(int in)
    {
        far_plane = in;
        dirty = true;
    }

    int get_raycast_threshold() const
    {
        return raycast_threshold;
    }

    void set_raycast_threshold(int in)
    {
        raycast_threshold = in;
        dirty = true;
    }

    void set_eye(glm::dvec3 in)
    {
        assert(is_real(in));
        dirty = true;
        eye = in;
    }

    void inc_eye(glm::dvec3 deye)
    {
        set_eye(eye + deye);
    }

    float get_theta() const
    {
        return theta;
    }

    void set_theta(float in)
    {
        assert(is_real(in));
        dirty = true;
        theta = in;
    }

    void inc_theta(float dtheta)
    {
        set_theta(theta + dtheta);
    }

    float get_phi() const
    {
        return phi;
    }

    void set_phi(float in)
    {
        assert(is_real(in));
        dirty = true;
        phi = glm::clamp(in, 0.01f, 3.14f);
    }

    void inc_phi(float dphi)
    {
        set_phi(phi + dphi);
    }

    float get_fovy_radians() const
    {
        return fovy_radians;
    }

    void set_fovy_radians(float in)
    {
        assert(is_real(in));
        assert(in > 0);
        dirty = true;
        fovy_radians = in;
    }

    void set_window_size(int x, int y)
    {
        dirty = true;
        window_x = x;
        window_y = y;
    }

    bool get_fog() const
    {
        return fog_enabled;
    }

    void set_fog(bool in)
    {
        fog_enabled = in;
    }

    // Move by the specified multiples of the normal right, up, and
    // forward vectors respectively.
    void frenet_move(float right, float up, float forward)
    {
        fix_dirty();
        inc_eye(glm::dvec3(
            right * right_vector
          + up * up_vector
          + forward * forward_normal_vector));
    }

    // Getters for derived linear algebra data.
    void get_eye(glm::ivec3* out_group, glm::vec3* out_residue = nullptr)
    {
        fix_dirty();
        if (out_group) *out_group = eye_group;
        if (out_residue) *out_residue = eye_residue;
    }

    glm::vec3 get_forward_normal()
    {
        fix_dirty();
        return forward_normal_vector;
    }

    glm::vec3 get_right_normal()
    {
        fix_dirty();
        return right_vector;
    }

    glm::vec3 get_up_normal()
    {
        fix_dirty();
        return up_vector;
    }

    glm::mat4 get_residue_view()
    {
        fix_dirty();
        return residue_view_matrix;
    }

    glm::mat4 get_projection()
    {
        fix_dirty();
        return projection_matrix;
    }

    glm::mat4 get_residue_vp()
    {
        fix_dirty();
        return residue_vp_matrix;
    }
};

} // end namespace
#endif /* !MYRICUBE_CAMERA_HH_ */
