// I'm using "app" to mean "some program that plays around with the
// voxel world". This is the API it is expected to follow.

#ifndef MYRICUBE_APP_HH_
#define MYRICUBE_APP_HH_

#include "myricube.hh"

#include "chunk.hh"
#include "window.hh"

namespace myricube {

// Name of the app.
extern const char app_name[];

// Called once at the start; also an opportunity to add key targets.
void app_init(VoxelWorld& world, Window& window);

// Called per frame.
void app_update(VoxelWorld& world);

} // end namespace
#endif /* !MYRICUBE_APP_HH_ */
