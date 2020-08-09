# Makefile for cross-compiling FOR Windows 64-bit ON Linux using mingw-w64.
# (Makes my life easier).
CC=x86_64-w64-mingw32-gcc
CXX=x86_64-w64-mingw32-g++
CPPFLAGS=-I ../glfw/include -I ../glad/include -I ../include -I ../glm
CXXFLAGS=-Wall -Wextra -O3 -g -std=c++17
CFLAGS=  -Wall -Wextra -O3 -g

default: all

-include cckiss/Makefile

glfw-cmake:
	cd glfw-build && cmake -DCMAKE_TOOLCHAIN_FILE=CMake/x86_64-w64-mingw32.cmake ../../glfw

glfw-build/src/libglfw3.a: glfw-cmake
	cd glfw-build && make

OBJS=cckiss/src/myricube.cc.o \
     cckiss/src/window.cc.o \
     cckiss/src/renderer.cc.o \
     cckiss/glad/src/glad.c.o \
     glfw-build/src/libglfw3.a \

all: cckiss/apps/random-walk.cc.o \
     myricube-random-walk.exe \
     $(OBJS)

# Empty line for backslash

# myricube-planet.exe: $(OBJS) cckiss/apps/marlo-planet.cc.o
# 	$(CXX) $(OBJS) cckiss/apps/marlo-planet.cc.o -ldl -lGL -lpthread -o myricube-planet.exe


myricube-%.exe: $(OBJS) cckiss/apps/%.cc.o
	$(CXX) $(OBJS) cckiss/apps/$*.cc.o -static -o myricube-$*.exe -lgdi32
