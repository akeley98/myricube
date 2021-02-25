#! /usr/bin/env python3
# File that prints out the build.ninja file to build myricube.
# ALL FILES LISTED MUST USE ONLY THE CHARACTERS [A-Z][a-z][0-9]_-./

import os

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
    "src/DummyVk.cc",
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

# Linker arguments
linkerflags = "glfw-build/src/libglfw3.a -ldl -lpthread"


        ## PREFER TO MAKE CHANGES ABOVE THIS LINE ##


# C and C++ compilers
cc = os.environ.get("CC", "cc")
cxx = os.environ.get("CXX", "c++")

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

# Invoke rules for compiling files, place stuff in build/ directory.
for src in c_files: print("build build/%s.o: cc %s" % (src, src))
for src in cxx_files: print("build build/%s.o: cxx %s" % (src, src))

# Link files together.
print("""
rule link
  command = %s $in -o $out %s
""" % (cxx, linkerflags))

objects = ["build/%s.o" % src for src in (c_files + cxx_files)]
print("build myricube-gl-bin: link %s" % ' '.join(objects))
