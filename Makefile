default: myricube-bin
	./myricube

CPPFLAGS=-I glfw/include -I ./glad/include -I include -I glm
CXXFLAGS=-Wall -Wextra -O3 -g -fPIC -std=c++17
CFLAGS=  -Wall -Wextra -O3 -g -fPIC

-include cckiss/Makefile
-include OBJS-list

glfw-cmake:
	cd glfw-build && cmake ../glfw

glfw-build/src/libglfw3.a: glfw-cmake
	cd glfw-build && $(MAKE)

myricube-windows:
	cd windows64 && $(MAKE)

all: myricube-bin myricube-windows

myricube-bin: $(OBJS)
	$(CXX) $(OBJS) -ldl -lpthread -o myricube-bin
