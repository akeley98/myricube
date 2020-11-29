default: myricube-bin
	./myricube

CPPFLAGS=-I glfw/include -I ./glad/include -I include -I glm -I vk-helpers
CFLAGS=  -Wall -Wextra -Wno-missing-field-initializers -O3 -g -fPIC
CXXFLAGS=$(CFLAGS) -std=c++17

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

myricube-bin:  $(OBJS) $(GL_OBJS) $(VK_OBJS) vk-shaders
	$(CXX) $(OBJS) $(GL_OBJS) $(VK_OBJS) -ldl -lpthread -lvulkan -o myricube-bin

myricube-vk-bin: $(OBJS) $(GL_DUMMY_OBJS) $(VK_OBJS) vk-shaders
	$(CXX)   $(OBJS) $(GL_DUMMY_OBJS) $(VK_OBJS) -ldl -lpthread -lvulkan -o myricube-vk-bin

myricube-gl-bin: $(OBJS) $(GL_OBJS) $(VK_DUMMY_OBJS)
	$(CXX)   $(OBJS) $(GL_OBJS) $(VK_DUMMY_OBJS) -ldl -lpthread -o myricube-gl-bin

# Rebuild vulkan shaders every time for now due to no method to track
# include dependencies.
vk-shaders:
	glslangValidator myricube-data/vk/mesh.vert -V -o myricube-data/vk/mesh.vert.spv
	glslangValidator myricube-data/vk/mesh.frag -V -o myricube-data/vk/mesh.frag.spv
	glslangValidator myricube-data/vk/raycast.vert -V -o myricube-data/vk/raycast.vert.spv
	glslangValidator myricube-data/vk/raycast.frag -V -o myricube-data/vk/raycast.frag.spv
	glslangValidator myricube-data/vk/background.vert -V -o myricube-data/vk/background.vert.spv
	glslangValidator myricube-data/vk/background.frag -V -o myricube-data/vk/background.frag.spv

py: libmyricube-cvoxel.so myricube-bin
	echo "NOTE: run './iMyricube.py [your-directory]' to change which world"
	python3 -i iMyricube.py myricube-data/PyWorld/

py-demo: libmyricube-cvoxel.so myricube-bin
	myricube_py_demo=1 python3 iMyricube.py myricube-data/PyDemoWorld/
