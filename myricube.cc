// Myricube: Experimental voxel renderer using hybrid raycast/mesh
// OpenGL rendering.
//
// David Zhao Akeley 2020 (not a great year)
//
// This was inspired by my discontent with Minecraft's low render
// distance, even on high settings.
//
// Like most voxel renderers, this renderer subdivides an "infinite"
// grid of voxels into cubic chunks. Nearby chunks are drawn using
// conventional methods (by converting the voxels into a mesh of
// triangles). Farther chunks are instead drawn using raycasting: a
// minimal AABB is calculated that contains all the solid voxels of
// the chunk, and the AABB itself is drawn using a raycasting fragment
// shader that checks for collisions with the voxels contained in the
// AABB. The voxel data itself (when raycasting is used) is stored
// using a 3D texture -- this is considerably more memory efficient
// than a triangle mesh.
//
// Cubes of chunks are organized into larger chunk groups. For OpenGL
// efficiency, chunks within the same chunk group share GPU resources
// (3D textures, buffer storage, etc.). Points in space are frequently
// expressed as "group" and "residue" coordinates -- these are the
// coordinates (in chunk-group-size units) of the lower-left of the
// chunk the point is in, and the remaining offset within the chunk
// (is that clear?). For example, if chunk groups were 100 x 100 x 100
// voxels, then point (1, 502.5, -1) has chunk coordinate (0, 5, -1)
// and residue (1, 2.5, 99).
//
// The primary benefits of raycasting are the aforementioned memory
// efficiency and the reduced number of triangle vertices processed
// (12 triangles are shared for all voxels per chunk, instead of up to
// 12 triangles per solid voxel). There are several drawbacks that make
// raycasting suitable only for distant chunks:
//
// * Raycasting cost is roughly linear in the size of the rendered
//   triangles (because the raycast fragment shader is the most
//   expensive part). This makes raycasting very expensive for large
//   nearby triangles.
//
// * Raycasting is more prone to rounding errors. This creates
//   unsightly gaps between voxels that are more obvious when nearby.
//
// * Raycasting causes "incorrect" writes to the depth buffer (because
//   the depth by default is that of the AABB, and not the actual
//   voxel), so raycast voxels may incorrectly occlude other geometry
//   (imagine a Minecraft mob standing on a voxel in the center of a
//   chunk). There are potential solutions but they might be expensive.
//
// * Circular reasoning, but in general I've designed the raycaster to
//   be cheap-but-ugly.
//
// This motivates the hybrid renderer I've come up with.

#include "myricube.hh"

#include <stdio.h>

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "window.hh"

namespace myricube {

// Absolute path of the executable, minus the -bin, plus -data/
// This is where shaders and stuff are stored.
std::string data_directory;

std::string expand_filename(const std::string& in)
{
    if (data_directory.size() == 0) {
        throw std::logic_error("Cannot call expand_filename before main");
    }
    return data_directory + in;
}

bool ends_with_dash_bin(const std::string& in)
{
    auto sz = in.size();
    return sz >= 4 and
           in[sz-4] == '-' and
           in[sz-3] == 'b' and
           in[sz-2] == 'i' and
           in[sz-1] == 'n';
           // sigh...
}

int Main(std::vector<std::string> args)
{
    if (args.at(0)[0] != '/') {
        fprintf(stderr, "%s should be absolute path\n"
            "(call through wrapper script).\n", args[0].c_str());
        return 1;
    }
    data_directory = args[0];
    // if (!data_directory.ends_with("-bin")) {
    if (!ends_with_dash_bin(data_directory)) {
        fprintf(stderr, "%s should end with '-bin'\n",
            args[0].c_str());
        return 1;
    }
    for (int i = 0; i < 4; ++i) data_directory.pop_back();
    data_directory += "-data/";

    Window window([] (int x, int y) { printf("%i %i\n", x, y); });
    while (window.update_swap_buffers(10)) continue;

    return 0;
}

} // end namespace

int main(int argc, char** argv)
{
    std::vector<std::string> args;
    for (int i = 0; i < argc; ++i) {
        args.emplace_back(argv[i]);
    }
    return myricube::Main(std::move(args));
}
