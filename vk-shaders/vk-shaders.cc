// C++ file containing the compiled SPIR-V shaders. The makefile
// contains the actual rule for converting GLSL to SPIR-V, then to the
// .xxd files containing comma-separated byte values.
// NOTE: For this to work, the makefile in this directory must be
// run and completed BEFORE cckiss compiles this C++ file.

#include <vector>
#include "vk-shaders.hh"

namespace myricube {

const std::vector<unsigned char> background_vert_spv {
#include "./background.vert.spv.xxd"
};

const std::vector<unsigned char> background_frag_spv {
#include "./background.frag.spv.xxd"
};

const std::vector<unsigned char> mesh_vert_spv {
#include "./mesh.vert.spv.xxd"
};

const std::vector<unsigned char> mesh_frag_spv {
#include "./mesh.frag.spv.xxd"
};

const std::vector<unsigned char> raycast_vert_spv {
#include "./raycast.vert.spv.xxd"
};

const std::vector<unsigned char> raycast_frag_spv {
#include "./raycast.frag.spv.xxd"
};

} // end namespace

