#include "window.hh"

#include <atomic>
#include <ctype.h>
#include <stdio.h>
#include <unordered_map>

#include "glad/glad.h"

static constexpr double fps_report_interval = 0.5;

namespace myricube {

// For now I'm just indicating mouse buttons with negative numbers.
// I'm using the X11 numbers I'm used to.
static inline int keycode_from_glfw_button(int button)
{
    switch (button) {
        case GLFW_MOUSE_BUTTON_LEFT: return -1;
        case GLFW_MOUSE_BUTTON_MIDDLE: return -2;
        case GLFW_MOUSE_BUTTON_RIGHT: return -3;
        case 3: return -8;
        case 4: return -9;
        case 5: return -10;
        case 6: return -11;
        default: return -12;
    }
}

// A pointer to the myricube Window object that created a glfw window
// is stored as the user pointer by the GLFW window. This function
// gets it back.
static inline Window& get_Window(GLFWwindow* window)
{
    Window* result = static_cast<Window*>(glfwGetWindowUserPointer(window));
    assert(result != nullptr);
    return *result;
}

Window::Window(OnWindowResize on_window_resize_)
{
    on_window_resize = std::move(on_window_resize_);

    // Should factor glfw initialization and hint setting out if I want
    // Window to be threadsafe (not sure if I really care for now...)
    if (!glfwInit()) {
        panic("Failed to initialize glfw");
    }
    glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    window = glfwCreateWindow(1920, 1080, "Myricube", nullptr, nullptr);

    if (window == nullptr) {
        panic("Could not initialize window");
    }

    glfwSetWindowUserPointer(window, this);
    glfwSetFramebufferSizeCallback(window, frame_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetCharModsCallback(window, character_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);

    gl_make_current();
    if (!gladLoadGL()) {
        panic("gladLoadGL failure");
    }
    glfwGetFramebufferSize(window, &frame_x, &frame_y);
    on_window_resize(frame_x, frame_y);
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

void Window::gl_make_current() const
{
    glfwMakeContextCurrent(window);
}

bool Window::update_events(float* out_dt)
{
    // Calculate dt. Require at least 1 ms between calls (workaround).
    double now;
    double dt;
    do {
        now = glfwGetTime();
        dt = now - previous_update;
    } while (dt < 0.001);

    previous_update = now;
    if (out_dt) *out_dt = dt;

    if (glfwWindowShouldClose(window)) return false;

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

    // Update FPS and frame time. NOTE: This is totally wrong now
    // that rendering is in a separate thread.
    ++frames;
    next_frame_time = std::max(next_frame_time, dt);
    if (now - previous_fps_update >= fps_report_interval) {
        fps = frames / (now - previous_fps_update);
        previous_fps_update = now;
        frames = 0;
        frame_time = next_frame_time;
        next_frame_time = 0;
    }
    return true;
}

void Window::swap_buffers()
{
    glfwSwapBuffers(window);
}

// Given that the key/mouse button with the specified key code has
// been pressed, search for a successful bind for it and call the
// KeyTarget for that bind. If successful, store the keycode and
// activated KeyTarget in the pressed_keys_map so we will do the right
// thing when the key is released (even if the key is re-bound in the
// meantime).
void Window::handle_down(int keycode, float amount)
{
    auto it = pressed_keys_map.find(keycode);
    KeyArg arg;
    arg.repeat = it != pressed_keys_map.end();
    arg.dt = 0.0;
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

void Window::frame_size_callback(GLFWwindow* window, int x, int y)
{
    Window& w = get_Window(window);
    if (w.on_window_resize) w.on_window_resize(x, y);
    // Invalidate cursor position on window resize.
    w.cursor_x = -1;
    w.cursor_y = -1;
}

void Window::key_callback(
    GLFWwindow* window, int key, int, int action, int)
{
    Window& w = get_Window(window);
    action != GLFW_RELEASE ? w.handle_down(key, 1) : w.handle_up(key);
}

void Window::character_callback(
    GLFWwindow*, unsigned int, int)
{
    // Might be useful later.
    // fprintf(stderr, "Codepoint input: %Xh   mods %i\n", codepoint, mods);
}

void Window::cursor_position_callback(
    GLFWwindow* window, double xpos, double ypos)
{
    Window& w = get_Window(window);
    bool valid = (w.cursor_x >= 0 and w.cursor_y >= 0);

    double dx = xpos - w.cursor_x;
    double dy = ypos - w.cursor_y;
    w.cursor_x = xpos;
    w.cursor_y = ypos;

    if (valid) {
        for (auto& pair : w.pressed_keys_map) {
            pair.second->mouse_rel_x += dx;
            pair.second->mouse_rel_y += dy;
        }
    }
}

void Window::mouse_button_callback(
    GLFWwindow* window, int button, int action, int)
{
    Window& w = get_Window(window);
    int keycode = keycode_from_glfw_button(button);
    action != GLFW_RELEASE ? w.handle_down(keycode, 1) : w.handle_up(keycode);
}

void Window::scroll_callback(GLFWwindow* window, double x, double y)
{
    Window& w = get_Window(window);
    // Again I'm just using the Unix mouse numbers I know.
    if (x < 0) {
        w.handle_down(-7, x);
        w.handle_up(-7);
    }
    if (x > 0) {
        w.handle_down(-6, x);
        w.handle_up(-6);
    }
    if (y < 0) {
        w.handle_down(-5, y);
        w.handle_up(-5);
    }
    if (y > 0) {
        w.handle_down(-4, y);
        w.handle_up(-4);
    }
}

static std::unordered_map<std::string, int> make_key_name_map();

int keycode_from_name(std::string name)
{
    static const std::unordered_map<std::string, int> m = make_key_name_map();
    for (char& c : name) {
        c = toupper(c);
        if (c == '-') c = '_';
    }
    auto it = m.find(name);
    return it != m.end() ? it->second : 0;
}

static std::unordered_map<std::string, int> make_key_name_map()
{
    std::unordered_map<std::string, int> m;

    m["MOUSE_1"] = -1;
    m["LEFT_MOUSE"] = -1;
    m["MOUSE_2"] = -2;
    m["MIDDLE_MOUSE"] = -2;
    m["MOUSE_3"] = -3;
    m["RIGHT_MOUSE"] = -3;
    m["MOUSE_4"] = -4;
    m["SCROLL_UP"] = -4;
    m["MOUSE_5"] = -5;
    m["SCROLL_DOWN"] = -5;
    m["MOUSE_6"] = -6;
    m["SCROLL_LEFT"] = -6;
    m["MOUSE_7"] = -7;
    m["SCROLL_RIGHT"] = -7;
    m["MOUSE_8"] = -8;
    m["THUMB_BUTTON"] = -8;
    m["X1"] = -8;
    m["MOUSE_9"] = -9;
    m["THUMB_BUTTON_2"] = -9;
    m["X2"] = -9;
    m["MOUSE_10"] = -10;
    m["MOUSE_11"] = -11;
    m["MOUSE_12"] = -12;

    m["SPACE"] = GLFW_KEY_SPACE;
    m["APOSTROPHE"] = GLFW_KEY_APOSTROPHE;
    m["COMMA"] = GLFW_KEY_COMMA;
    m["MINUS"] = GLFW_KEY_MINUS;
    m["PERIOD"] = GLFW_KEY_PERIOD;
    m["SLASH"] = GLFW_KEY_SLASH;
    m["0"] = GLFW_KEY_0;
    m["1"] = GLFW_KEY_1;
    m["2"] = GLFW_KEY_2;
    m["3"] = GLFW_KEY_3;
    m["4"] = GLFW_KEY_4;
    m["5"] = GLFW_KEY_5;
    m["6"] = GLFW_KEY_6;
    m["7"] = GLFW_KEY_7;
    m["8"] = GLFW_KEY_8;
    m["9"] = GLFW_KEY_9;
    m["SEMICOLON"] = GLFW_KEY_SEMICOLON;
    m["EQUAL"] = GLFW_KEY_EQUAL;
    m["A"] = GLFW_KEY_A;
    m["B"] = GLFW_KEY_B;
    m["C"] = GLFW_KEY_C;
    m["D"] = GLFW_KEY_D;
    m["E"] = GLFW_KEY_E;
    m["F"] = GLFW_KEY_F;
    m["G"] = GLFW_KEY_G;
    m["H"] = GLFW_KEY_H;
    m["I"] = GLFW_KEY_I;
    m["J"] = GLFW_KEY_J;
    m["K"] = GLFW_KEY_K;
    m["L"] = GLFW_KEY_L;
    m["M"] = GLFW_KEY_M;
    m["N"] = GLFW_KEY_N;
    m["O"] = GLFW_KEY_O;
    m["P"] = GLFW_KEY_P;
    m["Q"] = GLFW_KEY_Q;
    m["R"] = GLFW_KEY_R;
    m["S"] = GLFW_KEY_S;
    m["T"] = GLFW_KEY_T;
    m["U"] = GLFW_KEY_U;
    m["V"] = GLFW_KEY_V;
    m["W"] = GLFW_KEY_W;
    m["X"] = GLFW_KEY_X;
    m["Y"] = GLFW_KEY_Y;
    m["Z"] = GLFW_KEY_Z;
    m["LEFT_BRACKET"] = GLFW_KEY_LEFT_BRACKET;
    m["BACKSLASH"] = GLFW_KEY_BACKSLASH;
    m["RIGHT_BRACKET"] = GLFW_KEY_RIGHT_BRACKET;
    m["GRAVE_ACCENT"] = GLFW_KEY_GRAVE_ACCENT;
    m["WORLD_1"] = GLFW_KEY_WORLD_1;
    m["WORLD_2"] = GLFW_KEY_WORLD_2;

    m["ESCAPE"] = GLFW_KEY_ESCAPE;
    m["ENTER"] = GLFW_KEY_ENTER;
    m["TAB"] = GLFW_KEY_TAB;
    m["BACKSPACE"] = GLFW_KEY_BACKSPACE;
    m["INSERT"] = GLFW_KEY_INSERT;
    m["DELETE"] = GLFW_KEY_DELETE;
    m["RIGHT"] = GLFW_KEY_RIGHT;
    m["LEFT"] = GLFW_KEY_LEFT;
    m["DOWN"] = GLFW_KEY_DOWN;
    m["UP"] = GLFW_KEY_UP;
    m["PAGE_UP"] = GLFW_KEY_PAGE_UP;
    m["PAGE_DOWN"] = GLFW_KEY_PAGE_DOWN;
    m["HOME"] = GLFW_KEY_HOME;
    m["END"] = GLFW_KEY_END;
    m["CAPS_LOCK"] = GLFW_KEY_CAPS_LOCK;
    m["SCROLL_LOCK"] = GLFW_KEY_SCROLL_LOCK;
    m["NUM_LOCK"] = GLFW_KEY_NUM_LOCK;
    m["PRINT_SCREEN"] = GLFW_KEY_PRINT_SCREEN;
    m["PAUSE"] = GLFW_KEY_PAUSE;
    m["F1"] = GLFW_KEY_F1;
    m["F2"] = GLFW_KEY_F2;
    m["F3"] = GLFW_KEY_F3;
    m["F4"] = GLFW_KEY_F4;
    m["F5"] = GLFW_KEY_F5;
    m["F6"] = GLFW_KEY_F6;
    m["F7"] = GLFW_KEY_F7;
    m["F8"] = GLFW_KEY_F8;
    m["F9"] = GLFW_KEY_F9;
    m["F10"] = GLFW_KEY_F10;
    m["F11"] = GLFW_KEY_F11;
    m["F12"] = GLFW_KEY_F12;
    m["F13"] = GLFW_KEY_F13;
    m["F14"] = GLFW_KEY_F14;
    m["F15"] = GLFW_KEY_F15;
    m["F16"] = GLFW_KEY_F16;
    m["F17"] = GLFW_KEY_F17;
    m["F18"] = GLFW_KEY_F18;
    m["F19"] = GLFW_KEY_F19;
    m["F20"] = GLFW_KEY_F20;
    m["F21"] = GLFW_KEY_F21;
    m["F22"] = GLFW_KEY_F22;
    m["F23"] = GLFW_KEY_F23;
    m["F24"] = GLFW_KEY_F24;
    m["F25"] = GLFW_KEY_F25;
    m["KP_0"] = GLFW_KEY_KP_0;
    m["KP_1"] = GLFW_KEY_KP_1;
    m["KP_2"] = GLFW_KEY_KP_2;
    m["KP_3"] = GLFW_KEY_KP_3;
    m["KP_4"] = GLFW_KEY_KP_4;
    m["KP_5"] = GLFW_KEY_KP_5;
    m["KP_6"] = GLFW_KEY_KP_6;
    m["KP_7"] = GLFW_KEY_KP_7;
    m["KP_8"] = GLFW_KEY_KP_8;
    m["KP_9"] = GLFW_KEY_KP_9;
    m["KP_DECIMAL"] = GLFW_KEY_KP_DECIMAL;
    m["KP_DIVIDE"] = GLFW_KEY_KP_DIVIDE;
    m["KP_MULTIPLY"] = GLFW_KEY_KP_MULTIPLY;
    m["KP_SUBTRACT"] = GLFW_KEY_KP_SUBTRACT;
    m["KP_ADD"] = GLFW_KEY_KP_ADD;
    m["KP_ENTER"] = GLFW_KEY_KP_ENTER;
    m["KP_EQUAL"] = GLFW_KEY_KP_EQUAL;
    m["LEFT_SHIFT"] = GLFW_KEY_LEFT_SHIFT;
    m["LSHIFT"] = GLFW_KEY_LEFT_SHIFT;
    m["LEFT_CONTROL"] = GLFW_KEY_LEFT_CONTROL;
    m["LCTRL"] = GLFW_KEY_LEFT_CONTROL;
    m["LEFT_ALT"] = GLFW_KEY_LEFT_ALT;
    m["LALT"] = GLFW_KEY_LEFT_ALT;
    m["LEFT_SUPER"] = GLFW_KEY_LEFT_SUPER;
    m["LGUI"] = GLFW_KEY_LEFT_SUPER;
    m["RIGHT_SHIFT"] = GLFW_KEY_RIGHT_SHIFT;
    m["RSHIFT"] = GLFW_KEY_RIGHT_SHIFT;
    m["RIGHT_CONTROL"] = GLFW_KEY_RIGHT_CONTROL;
    m["RCTRL"] = GLFW_KEY_RIGHT_CONTROL;
    m["RIGHT_ALT"] = GLFW_KEY_RIGHT_ALT;
    m["RALT"] = GLFW_KEY_RIGHT_ALT;
    m["RIGHT_SUPER"] = GLFW_KEY_RIGHT_SUPER;
    m["RGUI"] = GLFW_KEY_RIGHT_SUPER;
    m["MENU"] = GLFW_KEY_MENU;

    return m;
}


} // end namespace
