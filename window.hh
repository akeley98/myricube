// Window class. Creates OpenGL context and manages user input via callbacks.
// Note: Although this is written as a class, I'm not confident this will
// work if there are multiple windows. Still, I avoid global state by habit.

#ifndef MYRICUBE_WINDOW_HH_
#define MYRICUBE_WINDOW_HH_

#include "myricube.hh"

#include <functional>
#include <unordered_map>
#include <vector>

#include "SDL2/SDL.h"

namespace myricube {

constexpr float max_dt = 1/15.f;

// Callback type of window resize handler (takes x,y of new window size).
using OnWindowResize = std::function<void(int,int)>;

struct KeyArg;

// Collection of callbacks that can be bound to a key. The KeyArg
// gives additional information.
struct KeyTarget
{
    // Called when the key is pressed. Return value is "success", if
    // not successful, another key binding may be tried. If this
    // function is empty, it is presumed successful.
    //
    // KeyArg::repeat indicates whether this is a repeat key.
    std::function<bool(KeyArg)> down;

    // Called when the key is released. Return code is meaningless for now.
    std::function<bool(KeyArg)> up;

    // Called every frame that the key is down, with KeyArg::dt
    // being the elapsed seconds of this frame. This is capped to
    // max_dt. Return code is meaningless for now.
    std::function<bool(KeyArg)> per_frame;

    // Ignore me.
    float mouse_rel_x = 0;
    float mouse_rel_y = 0;
};

struct KeyArg
{
    // Indicates key repeat if true.
    bool repeat = false;
    // When applicable, seconds since the previous frame (subject to
    // the max_dt cap).
    float dt = 0.0f;
    // Mystery amount, might be more useful later. Sort of intended to
    // indicate the "strength" of this operation.
    float amount = 1.0f;

    // Relative position of the mouse cursor compared to where it was
    // when this key was first pressed.
    float mouse_rel_x = 0;
    float mouse_rel_y = 0;
};

class Window
{
    // SDL window pointer and OpenGL context.
    SDL_Window* window = nullptr;
    void* context = nullptr;

    // Current size in pixels of the window.
    int window_x = 1280, window_y = 960;

    // Mapping from names of key targets to the actual key target structures.
    std::unordered_map<std::string, KeyTarget> key_target_map;

    // Mapping from key codes (or mouse buttons) to names of key
    // targets they can call (in reverse priority order). Only the
    // first successful keybinding is called.
    std::unordered_map<int, std::vector<std::string>> keycode_map;

    // A keycode is a key in this map iff it is pressed and it
    // resolved to a successful KeyTarget, which is pointed-to by the
    // value corresponding to the keycode key.
    std::unordered_map<int, KeyTarget*> pressed_keys_map;

    // Window resize callback.
    OnWindowResize on_window_resize;

    // milliseconds (since SDL initialization?) of previous
    // update_swap_buffers call.
    int64_t previous_update_ms = SDL_GetTicks();

  public:
    // Construct the window with a callback that is called when the
    // window is resized.
    Window(OnWindowResize);
    ~Window();
    Window(Window&&) = delete;

    // Set window title.
    void set_title(const std::string& title);

    // Get window size as pair of ints.
    void get_window_size(int* x, int* y) const
    {
        if (x) *x = window_x;
        if (y) *y = window_y;
    }

    // Make the OpenGL context of this window current.
    void gl_make_current();

    // Provide a name for the given key target. Physical keys can then
    // be bound to this name.
    void add_key_target(std::string target_name, KeyTarget key_target)
    {
        key_target_map.emplace(
            std::move(target_name),
            std::move(key_target));
    }

    // Bind this keycode to the key target with the given name.
    void bind_keycode(int keycode, std::string target_name)
    {
        keycode_map[keycode].emplace_back(std::move(target_name));
    }

    // Update events and swap OpenGL buffers. Return true iff the user
    // hasn't ordered the window closed yet. Parameter is the minimum
    // number of milliseconds between this call's completion and
    // the previous call's completion.
    bool update_swap_buffers(int64_t min_ms);

  private:
    void handle_down(int, float, float);
    void handle_up(int);
    void handle_key_down(const SDL_KeyboardEvent&, float);
    void handle_key_up(const SDL_KeyboardEvent&);
    void handle_mouse_wheel(const SDL_MouseWheelEvent&);
    void handle_mouse_down(const SDL_MouseButtonEvent&, float);
    void handle_mouse_up(const SDL_MouseButtonEvent&);
    void handle_mouse_motion(const SDL_MouseMotionEvent&);
    void handle_window_event(const SDL_WindowEvent&);
};

} // end namespace
#endif /* !MYRICUBE_WINDOW_HH_ */
