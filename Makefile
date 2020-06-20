default: run

CXX=clang++
CC=clang
CPPFLAGS=-I ./glad/include
CXXFLAGS=-Wall -Wextra -O3 -g -std=c++17
CFLAGS=  -Wall -Wextra -O3 -g

-include cckiss/Makefile

run: myricube-bin
	./myricube

OBJS=cckiss/myricube.cc.s \
     cckiss/window.cc.s \
     cckiss/renderer.cc.s \
     cckiss/glad/src/glad.c.s \
# Empty line for backslash

myricube-bin: $(OBJS)
	$(CXX) $(OBJS) -ldl -lGL -lSDL2 -o myricube-bin
