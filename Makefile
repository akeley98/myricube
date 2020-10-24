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

LIBOBJS=cckiss/src/cvoxel.cc.o cckiss/src/voxels.cc.o
libmyricube-cvoxel.so: $(LIBOBJS)
	$(CXX) $(LIBOBJS) -shared -o libmyricube-cvoxel.so

all: myricube-bin myricube-windows

myricube-bin: $(OBJS)
	$(CXX) $(OBJS) -ldl -lpthread -o myricube-bin

py: libmyricube-cvoxel.so myricube-bin
	echo "NOTE: run './iMyricube.py [your-directory]' to change which world"
	python3 -i iMyricube.py myricube-data/PyWorld/

py-demo: libmyricube-cvoxel.so myricube-bin
	myricube_py_demo=1 python3 iMyricube.py myricube-data/PyDemoWorld/
