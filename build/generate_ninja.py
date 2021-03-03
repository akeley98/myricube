#! /usr/bin/env python3
# File that prints out the build.ninja file to build myricube.
# ALL FILES LISTED MUST USE ONLY THE CHARACTERS [A-Z][a-z][0-9]_-./

import sys

# Include directories
include_dirs = (
    "include",
    "glfw/include",
    "glad/include",
    "glm",
    "vk-helpers",
    "vk-shaders",
)

# C files to compile and link into myricube
c_files = (
    "glad/src/glad.c",
)

# C++ files to compile and link into myricube
cxx_files = (
    "src/myricube.cc",
    "src/app.cc",
    "src/map.cc",
    "src/voxels.cc",
    "src/window.cc",
    "apps/AxisTest.cc",
    "apps/Congestion.cc",
    "apps/MarloPlanet.cc",
    "apps/MengerSponge.cc",
    "apps/RandomWalk.cc",
    "apps/Skygrid.cc",
    "apps/ViewWorld.cc",
    "apps/FastNoise.cpp",
    "src/RendererGL.cc",
    "src/RendererVk.cc",
    "vk-helpers/nvh/nvprint.cpp",
    "vk-helpers/nvvk/context_vk.cpp",
    "vk-helpers/nvvk/debug_util_vk.cpp",
    "vk-helpers/nvvk/error_vk.cpp",
    "vk-helpers/nvvk/extensions_vk.cpp",
    "vk-helpers/nvvk/swapchain_vk.cpp",
)

# GLSL files to compile and place in myricube-data
glsl_files = (
    "vk-shaders/mesh.vert",
    "vk-shaders/mesh.frag",
    "vk-shaders/mesh.mesh",
    "vk-shaders/mesh.task",
    "vk-shaders/raycast.vert",
    "vk-shaders/raycast.frag",
    "vk-shaders/background.vert",
    "vk-shaders/background.frag",
)

# Flags passed to C compiler
cflags = ' '.join((
    "-Wall",
    "-Wextra",
    "-Wno-missing-field-initializers",
    "-O3",
    "-g",
    "-fPIC",
))

# Flags passed to C++ compiler
cxxflags = cflags + " -std=c++17"

# Flags passed to GLSL compiler
glslflags = "-g"

# Linker arguments
linkerflags = "glfw-build/src/libglfw3.a -ldl -lpthread -lvulkan"


        ## PREFER TO MAKE CHANGES ABOVE THIS LINE ##


try:
    import userconfig
    config = userconfig.__dict__.get("config")
    if type(config) is not dict:
        raise TypeError("userconfig.py must provide a dict named config.")
except ImportError:
    print("userconfig.py not found, using defaults.", file=sys.stderr)
    config = {}

# C, C++, and glsl compilers, use defaults if not specified by user
cc = config.get("cc", "cc")
cxx = config.get("cxx", "c++")
glslc = config.get("glslc", "glslc")

# Additional flags provided by the user
cflags = ' '.join((cflags, config.get("cflags", "")))
cxxflags = ' '.join((cxxflags, config.get("cxxflags", "")))
glslflags = ' '.join((glslflags, config.get("glslflags", "")))

# Include path arguments
cppflags = "-I " + " -I ".join(include_dirs)

# Rule for C files
print("""
rule cc
  depfile = $out.d
  command = %s %s %s $in -c -o $out -MD -MF $out.d
""" % (cc, cppflags, cflags))

# Rule for C++ files
print("""
rule cxx
  depfile = $out.d
  command = %s %s %s $in -c -o $out -MD -MF $out.d
""" % (cxx, cppflags, cxxflags))

# Rule for glslc files
print("""
rule glsl
  depfile = $out.d
  command = %s %s %s $in -c -o $out -MD -MF $out.d
""" % (glslc, cppflags, glslflags))

# Invoke rules for compiling files, place C/C++ stuff in build/
# directory, GLSL stuff in myricube-data/
for src in c_files: print("build build/%s.o: cc %s" % (src, src))
for src in cxx_files: print("build build/%s.o: cxx %s" % (src, src))
for src in glsl_files: print("build windows64/myricube-data/%s.spv: glsl %s"
                             % (src, src))

# Link files together, explicit dependency on C/C++ object files
# and order-only dependencies on SPIR-V files.
print("""
rule link
  command = %s $in -o $out %s
""" % (cxx, linkerflags))

objects = ["build/%s.o" % src for src in (c_files + cxx_files)]
shaders = ["windows64/myricube-data/%s.spv" % src for src in glsl_files]
print("build myricube-bin: link %s || %s"
      % (' '.join(objects), ' '.join(shaders)))
