# Makefile for cross-compiling FOR Windows 64-bit ON Linux using
# mingw-w64. (Makes my life easier). This is pretty hacked-together
# to work on my own machine only.
CC=x86_64-w64-mingw32-gcc
CXX=x86_64-w64-mingw32-g++
CPPFLAGS=-I /usr/local/include -I ../glfw/include -I ../glad/include -I ../include -I ../glm -I ./vk-include -I ../vk-helpers -I vk-shaders
CXXFLAGS=-Wall -Wextra -O3 -g -std=c++17 -Wno-missing-field-initializers -municode
CFLAGS=  -Wall -Wextra -O3 -g -municode

VK_LIB=~/.wine/drive_c/VulkanSDK/1.2.154.1/Lib/vulkan-1.lib

default: myricube.exe myricube-gl.exe myricube-cvoxel.dll

-include cckiss/Makefile
-include ../OBJS-list

glfw-cmake:
	cd glfw-build && cmake -DCMAKE_TOOLCHAIN_FILE=CMake/x86_64-w64-mingw32.cmake ../../glfw

glfw-build/src/libglfw3.a: glfw-cmake
	cd glfw-build && $(MAKE)

W64_OPTIONS=-municode -static -lgdi32

# Note: I embed compiled SPIR-V shaders into the .exe file. Vulkan
# builds depend on compile-vk-shaders from parent directory makefile
# to be run. (see myricube-windows target from said makefile).
# There might be a race condition if I try to compile glsl here too, unsure.

myricube.exe:  $(OBJS) $(GL_OBJS) $(VK_OBJS)
	$(CXX) $(OBJS) $(GL_OBJS) $(VK_OBJS) $(VK_LIB) $(W64_OPTIONS) -o myricube.exe

myricube-vk.exe: $(OBJS) $(GL_DUMMY_OBJS) $(VK_OBJS)
	$(CXX)   $(OBJS) $(GL_DUMMY_OBJS) $(VK_OBJS) $(VK_LIB) $(W64_OPTIONS) -o myricube-vk.exe

myricube-gl.exe: $(OBJS) $(GL_OBJS) $(VK_DUMMY_OBJS)
	$(CXX)   $(OBJS) $(GL_OBJS) $(VK_DUMMY_OBJS) $(W64_OPTIONS) -o myricube-gl.exe

myricube-cvoxel.dll: $(LIBOBJS)
	$(CXX) $(LIBOBJS) -municode -static -shared -o myricube-cvoxel.dll
