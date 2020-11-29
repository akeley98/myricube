// Just view the world given by the myricube_world environment variable.
#include "app.hh"

#include <stdlib.h>

#include "EnvVar.hh"

namespace myricube {

class ViewWorld : public App
{
    VoxelWorld world { get_world_filename() };

    static const filename_string get_world_filename()
    {
        EnvVarFilename env("myricube_world", "");
        filename_string result = env;
        if (result.empty()) {
            throw std::runtime_error(
                "Empty or missing myricube_world environment variable");
        }
        return result;
    }

  public:
    VoxelWorld& update(float) override
    {
        return world;
    }
};

MYRICUBE_ADD_APP(ViewWorld)

} // end namespace
