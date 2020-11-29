# Makefile for cross-compiling FOR Windows 64-bit ON Linux using
# mingw-w64. (Makes my life easier). This is pretty hacked-together
# to work on my own machine only.
CC=x86_64-w64-mingw32-gcc
CXX=x86_64-w64-mingw32-g++
CPPFLAGS=-I ../glfw/include -I ../glad/include -I ../include -I ../glm -I ./vk-include -I ../vk-helpers
CXXFLAGS=-Wall -Wextra -O3 -g -std=c++17 -Wno-missing-field-initializers -municode
CFLAGS=  -Wall -Wextra -O3 -g -municode

VK_LIB=~/.wine/drive_c/VulkanSDK/1.2.154.1/Lib/vulkan-1.lib

default: myricube.exe myricube-gl.exe

-include cckiss/Makefile
-include ../OBJS-list

glfw-cmake:
	cd glfw-build && cmake -DCMAKE_TOOLCHAIN_FILE=CMake/x86_64-w64-mingw32.cmake ../../glfw

glfw-build/src/libglfw3.a: glfw-cmake
	cd glfw-build && $(MAKE)

W64_OPTIONS=-municode -static -lgdi32

myricube.exe:  $(OBJS) $(GL_OBJS) $(VK_OBJS) vk-shaders
	$(CXX) $(OBJS) $(GL_OBJS) $(VK_OBJS) $(VK_LIB) $(W64_OPTIONS) -o myricube.exe

myricube-vk.exe: $(OBJS) $(GL_DUMMY_OBJS) $(VK_OBJS) vk-shaders
	$(CXX)   $(OBJS) $(GL_DUMMY_OBJS) $(VK_OBJS) $(VK_LIB) $(W64_OPTIONS) -o myricube-vk.exe

myricube-gl.exe: $(OBJS) $(GL_OBJS) $(VK_DUMMY_OBJS)
	$(CXX)   $(OBJS) $(GL_OBJS) $(VK_DUMMY_OBJS) $(W64_OPTIONS) -o myricube-gl.exe


# Rebuild vulkan shaders every time for now due to no method to track
# include dependencies. TODO share this with other Makefile?
vk-shaders:
	glslangValidator myricube-data/vk/mesh.vert -V -o myricube-data/vk/mesh.vert.spv
	glslangValidator myricube-data/vk/mesh.frag -V -o myricube-data/vk/mesh.frag.spv
	glslangValidator myricube-data/vk/raycast.vert -V -o myricube-data/vk/raycast.vert.spv
	glslangValidator myricube-data/vk/raycast.frag -V -o myricube-data/vk/raycast.frag.spv
	glslangValidator myricube-data/vk/background.vert -V -o myricube-data/vk/background.vert.spv
	glslangValidator myricube-data/vk/background.frag -V -o myricube-data/vk/background.frag.spv
