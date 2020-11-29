// Just view the world given by the myricube_world environment variable.
#include "app.hh"

#include <stdlib.h>

namespace myricube {

class ViewWorld : public App
{
    VoxelWorld world { get_world_filename() };

#ifdef MYRICUBE_WINDOWS
    static const wchar_t* get_world_filename()
    {
        const wchar_t* result = _wgetenv(L"myricube_world");
        if (!result) {
            panic("Missing myricube_world environment variable.");
        }
        return result;
    }
#else
    static const char* get_world_filename()
    {
        const char* result = getenv("myricube_world");
        if (!result) {
            panic("Missing myricube_world environment variable.");
        }
        return result;
    }
#endif

  public:
    VoxelWorld& update(float) override
    {
        return world;
    }
};

MYRICUBE_ADD_APP(ViewWorld)

} // end namespace
