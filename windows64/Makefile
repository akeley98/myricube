# Makefile for cross-compiling FOR Windows 64-bit ON Linux using mingw-w64.
# (Makes my life easier).
CC=x86_64-w64-mingw32-gcc
CXX=x86_64-w64-mingw32-g++
CPPFLAGS=-I ../glfw/include -I ../glad/include -I ../include -I ../glm
CXXFLAGS=-Wall -Wextra -O3 -g -std=c++17
CFLAGS=  -Wall -Wextra -O3 -g

default: myricube.exe

-include cckiss/Makefile
-include ../OBJS-list

glfw-cmake:
	cd glfw-build && cmake -DCMAKE_TOOLCHAIN_FILE=CMake/x86_64-w64-mingw32.cmake ../../glfw

glfw-build/src/libglfw3.a: glfw-cmake
	cd glfw-build && $(MAKE)

myricube.exe: $(OBJS)
	$(CXX) $(OBJS) -static -o myricube.exe -lgdi32
