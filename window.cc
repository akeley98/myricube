#include "window.hh"

#include <stdio.h>

#include "glad/glad.h"
#include "SDL2/SDL.h"
#include "SDL2/SDL_opengl.h"

namespace myricube {

// For now I'm just indicating mouse buttons with negative numbers.
// I'm using the X11 numbers I'm used to.
int keycode_from_button_event(const SDL_MouseButtonEvent& e)
{
    switch (e.button) {
        default: return -10;
        case SDL_BUTTON_LEFT: return -1;
        case SDL_BUTTON_MIDDLE: return -2;
        case SDL_BUTTON_RIGHT: return -3;
        case SDL_BUTTON_X1: return -8;
        case SDL_BUTTON_X2: return -9;
    }
}

Window::Window(OnWindowResize on_window_resize_)
{
    on_window_resize = on_window_resize_;
    window = SDL_CreateWindow(
        "",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        window_x, window_y,
        SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE
    );
    if (window == nullptr) {
        panic("Could not initialize window", SDL_GetError());
    }
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 5);

    context = SDL_GL_CreateContext(window);
    if (!context) {
        panic("OpenGL context creation error:", SDL_GetError());
    }
    gl_make_current();

    if (!gladLoadGL()) {
        panic("gladLoadGL failure", "gladLoadGL failure");
    }
}

Window::~Window()
{
    fprintf(stderr, "Maybe I should write a destructor for Window...\n");
}

void Window::set_title(const std::string& title)
{
    SDL_SetWindowTitle(window, title.c_str());
}

void Window::gl_make_current()
{
    if (SDL_GL_MakeCurrent(window, context) < 0) {
        panic("SDL OpenGL context error", SDL_GetError());
    }
}

bool Window::update_swap_buffers(int64_t min_ms)
{
    // Calculate dt & wait for min_ms elapsed before previous swap.
    int64_t current_tick = 0;
    int64_t dt_ms = 0;
    for (; dt_ms < min_ms; SDL_Delay(1)) {
        current_tick = SDL_GetTicks();
        dt_ms = current_tick - previous_update_ms;
    }
    if (dt_ms > next_frame_time) next_frame_time = dt_ms;
    previous_update_ms += dt_ms;
    float dt = 1/1000.f * dt_ms;

    // Swap buffers as promised.
    SDL_GL_SwapWindow(window);

    // Call per-frame callbacks of pressed keys.
    KeyArg arg;
    arg.dt = std::min(dt, max_dt);
    for (auto& pair : pressed_keys_map) {
        arg.mouse_rel_x = pair.second->mouse_rel_x;
        arg.mouse_rel_y = pair.second->mouse_rel_y;
        auto& cb = pair.second->per_frame;
        if (cb) cb(arg);
    }

    // Check for new events.
    bool quit = false;
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        switch (event.type) {
            case SDL_KEYDOWN: handle_key_down(event.key, dt); break;
            case SDL_KEYUP: handle_key_up(event.key); break;
            case SDL_MOUSEWHEEL: handle_mouse_wheel(event.wheel); break;
            case SDL_MOUSEBUTTONDOWN: handle_mouse_down(event.button, dt); break;
            case SDL_MOUSEBUTTONUP: handle_mouse_up(event.button); break;
            case SDL_MOUSEMOTION: handle_mouse_motion(event.motion); break;
            case SDL_WINDOWEVENT: handle_window_event(event.window); break;
            case SDL_QUIT: quit = true; break;
        }
    }

    // Update FPS and frame time.
    ++frames;
    if (current_tick - previous_fps_update >= fps_interval_ms) {
        fps = 1000.0 * frames
            / (current_tick - previous_fps_update);
        previous_fps_update = current_tick;
        frames = 0;
        frame_time_ms = next_frame_time;
        next_frame_time = 0;
    }
    return !quit;
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

void Window::handle_key_down(const SDL_KeyboardEvent& e, float dt)
{
    handle_down(e.keysym.scancode, dt, 1.0f);
}

void Window::handle_key_up(const SDL_KeyboardEvent& e)
{
    handle_up(e.keysym.scancode);
}

void Window::handle_mouse_down(const SDL_MouseButtonEvent& e, float dt)
{
    handle_down(keycode_from_button_event(e), dt, 1.0f);
}

void Window::handle_mouse_up(const SDL_MouseButtonEvent& e)
{
    handle_up(keycode_from_button_event(e));
}

void Window::handle_mouse_wheel(const SDL_MouseWheelEvent& e)
{
    // Again I'm just using the Unix mouse numbers I know.
    // x seems to be backwards from what I expect in SDL
    if (e.x < 0) {
        handle_down(-7, 0, e.x);
        handle_up(-7);
    }
    if (e.x > 0) {
        handle_down(-6, 0, e.x);
        handle_up(-6);
    }
    if (e.y < 0) {
        handle_down(-5, 0, e.y);
        handle_up(-5);
    }
    if (e.y > 0) {
        handle_down(-4, 0, e.y);
        handle_up(-4);
    }
}

void Window::handle_mouse_motion(const SDL_MouseMotionEvent& e)
{
    float dx = e.xrel;
    float dy = e.yrel;
    for (auto& pair : pressed_keys_map) {
        pair.second->mouse_rel_x += dx;
        pair.second->mouse_rel_y += dy;
    }
}

void Window::handle_window_event(const SDL_WindowEvent& e)
{
    if (e.event == SDL_WINDOWEVENT_SIZE_CHANGED ||
        e.event == SDL_WINDOWEVENT_RESIZED) {
            window_x = e.data1;
            window_y = e.data2;
    }
    assert(on_window_resize);
    on_window_resize(window_x, window_y);
}

} // end namespace
