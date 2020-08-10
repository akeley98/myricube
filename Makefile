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

# Target all-w64 to REALLY make everything for Windows and Linux.
all-w64: all
	cd windows64 && $(MAKE)

myricube-bin: $(OBJS)
	$(CXX) $(OBJS) -ldl -lGL -lpthread -o myricube-bin
