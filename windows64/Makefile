# Makefile for cross-compiling FOR Windows 64-bit ON Linux using
# mingw-w64. (Makes my life easier). This is pretty hacked-together
# to work on my own machine only.
CC=x86_64-w64-mingw32-gcc
CXX=x86_64-w64-mingw32-g++
CPPFLAGS=-I ../glfw/include -I ../glad/include -I ../include -I ../glm -I ./vk-include -I ../vk-helpers
CXXFLAGS=-Wall -Wextra -O3 -g -std=c++17 -Wno-missing-field-initializers -municode
CFLAGS=  -Wall -Wextra -O3 -g -municode

VK_LIB=~/.wine/drive_c/VulkanSDK/1.2.154.1/Lib/vulkan-1.lib

default: myricube.exe

-include cckiss/Makefile
-include ../OBJS-list

glfw-cmake:
	cd glfw-build && cmake -DCMAKE_TOOLCHAIN_FILE=CMake/x86_64-w64-mingw32.cmake ../../glfw

glfw-build/src/libglfw3.a: glfw-cmake
	cd glfw-build && $(MAKE)

myricube.exe: myricube-gl-vk

W64_OPTIONS=-municode -static -o myricube.exe -lgdi32

myricube-gl-vk: $(OBJS) $(GL_OBJS) $(VK_OBJS) vk-shaders
	$(CXX)  $(OBJS) $(GL_OBJS) $(VK_OBJS) $(VK_LIB) $(W64_OPTIONS) -lvulkan

myricube-vk:   $(OBJS) $(GL_DUMMY_OBJS) $(VK_OBJS) vk-shaders
	$(CXX) $(OBJS) $(GL_DUMMY_OBJS) $(VK_OBJS) $(VK_LIB) $(W64_OPTIONS) -lvulkan

myricube-gl:   $(OBJS) $(GL_OBJS) $(VK_DUMMY_OBJS)
	$(CXX) $(OBJS) $(GL_OBJS) $(VK_DUMMY_OBJS) $(W64_OPTIONS)
