#include "window.hh"

#include <atomic>
#include <stdio.h>

#include "glad/glad.h"

namespace myricube {

// For now I'm just indicating mouse buttons with negative numbers.
// I'm using the X11 numbers I'm used to.
// int keycode_from_button_event(const SDL_MouseButtonEvent& e)
// {
//     switch (e.button) {
//         default: return -10;
//         case SDL_BUTTON_LEFT: return -1;
//         case SDL_BUTTON_MIDDLE: return -2;
//         case SDL_BUTTON_RIGHT: return -3;
//         case SDL_BUTTON_X1: return -8;
//         case SDL_BUTTON_X2: return -9;
//     }
// }

Window::Window(OnWindowResize on_window_resize_)
{
    on_window_resize = on_window_resize_;

    // Should factor glfw initialization and hint setting out if I want
    // Window to be threadsafe (not sure if I really care for now...)
    if (!glfwInit()) {
        panic("Failed to initialize glfw");
    }
    glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    window = glfwCreateWindow(window_x, window_y, "", nullptr, nullptr);

    if (window == nullptr) {
        panic("Could not initialize window");
    }
    glfwSetWindowUserPointer(window, this);
    gl_make_current();

    if (!gladLoadGL()) {
        panic("gladLoadGL failure");
    }
    on_window_resize(window_x, window_y);
}

Window::~Window()
{
    glfwDestroyWindow(window);
    window = nullptr;
}

void Window::set_title(const std::string& title)
{
    glfwSetWindowTitle(window, title.c_str());
}

void Window::set_title(const char* title)
{
    glfwSetWindowTitle(window, title);
}

void Window::gl_make_current()
{
    glfwMakeContextCurrent(window);
}

bool Window::update_swap_buffers()
{
    if (glfwWindowShouldClose(window)) return false;

    // Calculate dt.
    double now = glfwGetTime();
    double dt = now - previous_update;
    previous_update = now;

    // Swap buffers as promised.
    glfwSwapBuffers(window);

    // Call per-frame callbacks of pressed keys.
    glfwPollEvents();
    KeyArg arg;
    arg.dt = std::min(float(dt), float(max_dt));
    for (auto& pair : pressed_keys_map) {
        arg.mouse_rel_x = pair.second->mouse_rel_x;
        arg.mouse_rel_y = pair.second->mouse_rel_y;
        auto& cb = pair.second->per_frame;
        if (cb) cb(arg);
    }

    // Update FPS and frame time.
    ++frames;
    frame_time = std::max(frame_time, dt);
    if (now - previous_fps_update >= fps_interval) {
        fps = frames / (now - previous_fps_update);
        previous_fps_update = now;
        frames = 0;
        frame_time = next_frame_time;
        next_frame_time = 0;
    }
    return true;
}

// Given that the key/mouse button with the specified key code has
// been pressed, search for a successful bind for it and call the
// KeyTarget for that bind. If successful, store the keycode and
// activated KeyTarget in the pressed_keys_map so we will do the right
// thing when the key is released (even if the key is re-bound in the
// meantime).
void Window::handle_down(int keycode, float dt, float amount)
{
    auto it = pressed_keys_map.find(keycode);
    KeyArg arg;
    arg.repeat = it != pressed_keys_map.end();
    arg.dt = dt;
    arg.amount = amount;

    // In the repeat case, just recycle the cached KeyTarget ptr in the
    // pressed_keys map.
    if (arg.repeat) {
        arg.mouse_rel_x = it->second->mouse_rel_x;
        arg.mouse_rel_y = it->second->mouse_rel_y;
        auto& on_down = it->second->down;
        if (on_down) on_down(arg);
        return;
    }
    // Otherwise, we have to try to find a bind (string name) and convert
    // it to a KeyTarget to call.

    // Are there any binds for this physical key?
    auto keycode_targets_it = keycode_map.find(keycode);
    if (keycode_targets_it == keycode_map.end()) return;

    // If there are, try all the binds in reverse (priority) order,
    // remembering in the target ptr the first successful
    // one. (nullptr if there are no successful ones).
    const std::vector<std::string>& bind_strs = keycode_targets_it->second;
    KeyTarget* target = nullptr;
    for (auto str_it = rbegin(bind_strs); str_it != rend(bind_strs); ++str_it)
    {
        auto kt_it = key_target_map.find(*str_it);
        if (kt_it == key_target_map.end()) {
            fprintf(stderr, "No KeyTarget for '%s'\n", str_it->c_str());
            continue;
        }
        target = &(kt_it->second);
        if (! target->down) break;
        bool okay = kt_it->second.down(arg);
        if (okay) break;
        else {
            // Try next bind on failure.
            fprintf(stderr, "KeyTarget for '%s' failed\n", str_it->c_str());
            target = nullptr;
        }
    }

    // As per the comment about the pressed_keys_map, only if there is
    // a successful bind, remember it in the pressed_keys
    // dictionary. Also, we have to reset the mouse_rel_x/y counter,
    // which is slimily hidden in the KeyTarget as well.
    if (target) {
        target->mouse_rel_x = 0;
        target->mouse_rel_y = 0;
        pressed_keys_map[keycode] = target;
    }
}

// Look in the pressed keys map for the correct "on key up" callback to
// call, then remove the entry from the pressed keys map.
void Window::handle_up(int keycode)
{
    auto it = pressed_keys_map.find(keycode);
    if (it == pressed_keys_map.end()) {
        fprintf(stderr, "%i: no pressed_keys_map entry.\n", keycode);
        return;
    }

    KeyArg arg;
    arg.mouse_rel_x = it->second->mouse_rel_x;
    arg.mouse_rel_y = it->second->mouse_rel_y;
    if (it->second->up) {
        it->second->up(arg);
    }
    pressed_keys_map.erase(it);
}

} // end namespace
