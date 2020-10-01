// Camera class. Mostly just some vectors and matrices and stuff.

#ifndef MYRICUBE_CAMERA_HH_
#define MYRICUBE_CAMERA_HH_

#include "myricube.hh"

#include <assert.h>
#include <math.h>
#include <mutex>
#include <stddef.h>

namespace myricube {

// Matrices, etc. derived from a camera.
struct CameraTransforms
{
    // Eye coordinates.
    glm::dvec3 eye;

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

    // Near and far plane used to make the projection matrix.
    float near_plane;
    int far_plane;

    // Fog flags
    bool use_fog;
    bool use_black_fog;

    // Maximum number of new chunk groups added to GPU memory per frame.
    int max_frame_new_chunk_groups;

    int target_fragments;
    int screen_x, screen_y;
};

// Your typical camera info (eye position, far plane, fov, etc.).
//
// Used to communicate between render and control threads, so every
// member function is mutex protected.
class SyncCamera
{
    // Mutex, to be locked when reading/writing from this Camera.
    mutable std::mutex camera_mutex;

    // Position of the eye.
    glm::dvec3 eye = glm::dvec3(0.0, 0.0, 0.0);

    // Near and far plane for z-depth.
    // Far plane influences the chunk render distance.
    float near_plane = 0.1f;
    int far_plane = 512;

    // Horizontal and vertical angle camera is pointed in.
    float theta = 1.5707f, phi = 1.5707f;

    // Field of view (y direction), radians.
    float fovy_radians = 1.0f;

    // Window size in pixels.
    int window_x = 1, window_y = 1;

    // Maximum number of new chunk groups added to GPU memory per frame.
    int max_frame_new_chunk_groups = 10;

    // Fog setting.
    bool fog_enabled = true;
    bool black_fog = false;

    // Maximum number of fragments for the screen (integer
    // downsampling is done to meet this limit). Non-positive value
    // indicates no limit.
    int target_fragments = 0;

    // Frenet frame of the camera. This is updated every time the
    // camera angle changes (I don't try to be clever).
    glm::vec3 forward_normal_vector;
    glm::vec3 right_vector;
    glm::vec3 up_vector;

    void update_frenet_frame()
    {
        forward_normal_vector = glm::vec3(
            sinf(phi) * cosf(theta),
            cosf(phi),
            sinf(phi) * sinf(theta));

        right_vector = glm::normalize(
            glm::cross(forward_normal_vector, glm::vec3(0, 1, 0)));
        up_vector = glm::cross(right_vector, forward_normal_vector);
    }

  public:
    SyncCamera()
    {
        update_frenet_frame();
    }

    // Getters and setters for user-specified camera data.
    glm::dvec3 get_eye() const
    {
        std::lock_guard guard(camera_mutex);
        return eye;
    }

    float get_near_plane() const
    {
        std::lock_guard guard(camera_mutex);
        return near_plane;
    }

    void set_near_plane(float in)
    {
        std::lock_guard guard(camera_mutex);
        near_plane = in;
    }

    int get_far_plane() const
    {
        std::lock_guard guard(camera_mutex);
        return far_plane;
    }

    void set_far_plane(int in)
    {
        std::lock_guard guard(camera_mutex);
        far_plane = in;
    }

    void set_eye(glm::dvec3 in)
    {
        std::lock_guard guard(camera_mutex);
        assert(is_real(in));
        eye = in;
    }

    void inc_eye(glm::dvec3 deye)
    {
        set_eye(eye + deye);
    }

    float get_theta() const
    {
        std::lock_guard guard(camera_mutex);
        return theta;
    }

    void set_theta(float in)
    {
        std::lock_guard guard(camera_mutex);
        assert(is_real(in));
        theta = in;
        update_frenet_frame();
    }

    void inc_theta(float dtheta)
    {
        set_theta(theta + dtheta);
    }

    float get_phi() const
    {
        std::lock_guard guard(camera_mutex);
        return phi;
    }

    void set_phi(float in)
    {
        std::lock_guard guard(camera_mutex);
        assert(is_real(in));
        phi = glm::clamp(in, 0.01f, 3.14f);
        update_frenet_frame();
    }

    void inc_phi(float dphi)
    {
        set_phi(phi + dphi);
    }

    float get_fovy_radians() const
    {
        std::lock_guard guard(camera_mutex);
        return fovy_radians;
    }

    void set_fovy_radians(float in)
    {
        std::lock_guard guard(camera_mutex);
        assert(is_real(in));
        assert(in > 0);
        fovy_radians = in;
    }

    void set_window_size(int x, int y)
    {
        std::lock_guard guard(camera_mutex);
        window_x = x;
        window_y = y;
    }

    void get_window_size(int* x=nullptr, int* y=nullptr)
    {
        std::lock_guard guard(camera_mutex);
        if (x) *x = window_x;
        if (y) *y = window_y;
    }

    bool get_fog() const
    {
        std::lock_guard guard(camera_mutex);
        return fog_enabled;
    }

    void set_fog(bool in)
    {
        std::lock_guard guard(camera_mutex);
        fog_enabled = in;
    }

    bool use_black_fog() const
    {
        std::lock_guard guard(camera_mutex);
        return black_fog;
    }

    bool use_black_fog(bool in)
    {
        std::lock_guard guard(camera_mutex);
        return black_fog = in;
    }

    int get_max_frame_new_chunk_groups() const
    {
        std::lock_guard guard(camera_mutex);
        return max_frame_new_chunk_groups;
    }

    void set_max_frame_new_chunk_groups(int in)
    {
        std::lock_guard guard(camera_mutex);
        max_frame_new_chunk_groups = in;
    }

    int get_target_fragments() const
    {
        std::lock_guard guard(camera_mutex);
        return target_fragments;
    }

    void set_target_fragments(int in)
    {
        std::lock_guard guard(camera_mutex);
        target_fragments = in;
    }

    // Move by the specified multiples of the normal right, up, and
    // forward vectors respectively.
    void frenet_move(float right, float up, float forward)
    {
        std::lock_guard guard(camera_mutex);
        eye += glm::dvec3(
                    right * right_vector
                  + up * up_vector
                  + forward * forward_normal_vector);
    }

    // Convert camera data to actual camera transforms needed.
    explicit operator CameraTransforms() const
    {
        std::lock_guard guard(camera_mutex);
        CameraTransforms t;

        t.eye = eye;
        split_coordinate(eye, &t.eye_group, &t.eye_residue);
        t.forward_normal_vector = forward_normal_vector;
        t.right_vector = right_vector;
        t.up_vector = up_vector;

        t.residue_view_matrix = glm::lookAt(t.eye_residue,
                                            t.eye_residue + forward_normal_vector,
                                            glm::vec3(0, 1, 0));

        t.projection_matrix = glm::perspective(
            float(fovy_radians),
            float(window_x) / window_y,
            float(near_plane),
            float(far_plane));

        t.residue_vp_matrix = t.projection_matrix * t.residue_view_matrix;

        t.near_plane = near_plane;
        t.far_plane = far_plane;
        t.use_fog = fog_enabled;
        t.use_black_fog = black_fog;
        t.max_frame_new_chunk_groups = max_frame_new_chunk_groups;
        t.target_fragments = target_fragments;
        t.screen_x = window_x;
        t.screen_y = window_y;

        return t;
    }
};

} // end namespace
#endif /* !MYRICUBE_CAMERA_HH_ */
