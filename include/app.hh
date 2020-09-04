// I'm using "app" to mean "some program that plays around with the
// voxel world". Each is given some unique name and can be chosen and
// instantiated at runtime by name. This is the API apps are expected
// to follow.

#ifndef MYRICUBE_APP_HH_
#define MYRICUBE_APP_HH_

#include "myricube.hh"

#include "chunk.hh"
#include "window.hh"

namespace myricube {

// Apps should derive from this base class, providing a constructor
// that takes no arguments. If anything is to actually be drawn, a
// VoxelWorld needs to be instantiated and modified by the derived
// class (see chunk.hh -- not the best name in hindsight).
class App
{
  public:
    // Called every frame with (in theory) the frame time.  Returns a
    // reference to the current world to draw (this means the App is
    // responsible for the lifetime of the worlds it creates.)
    virtual VoxelWorld& update(float dt) = 0;

    // Optionally, override this function to provide more KeyTargets
    // that keys may be bound to.
    virtual void add_key_targets(Window&) {};

    // Apps may assume they are not moved through memory.
    App(App&&) = delete;
    App() = default;
    virtual ~App() = default;
};

// Then use this macro to add your app class by name to the global
// list of app factory functions.
#define MYRICUBE_ADD_APP(AppType) static AppNamer APP_NAMER##AppType( \
    #AppType, \
    __FILE__, \
    [] () -> App* { return new AppType(); } );

// Helper for above macro.
struct AppNamer
{
    AppNamer(std::string name, const char* filename, App* (*factory)());
};

// Run the factory function for the named app, or return null if no
// such app was linked-in.
App* new_named_app(const std::string& name);

void stderr_dump_app_names();

} // end namespace

#endif /* !MYRICUBE_APP_HH_ */
