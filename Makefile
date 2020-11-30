default: myricube-bin
	./myricube

CPPFLAGS=-I glfw/include -I ./glad/include -I include -I glm -I vk-helpers -I vk-shaders
CFLAGS=  -Wall -Wextra -Wno-missing-field-initializers -O3 -g -fPIC
CXXFLAGS=$(CFLAGS) -std=c++17

-include cckiss/Makefile
-include OBJS-list

glfw-cmake:
	cd glfw-build && cmake ../glfw

glfw-build/src/libglfw3.a: glfw-cmake
	cd glfw-build && $(MAKE)

libmyricube-cvoxel.so: $(LIBOBJS)
	$(CXX) $(LIBOBJS) -shared -o libmyricube-cvoxel.so

all: myricube-bin myricube-windows

# A bit of a hack: to make distribution easier, I generate a C++
# source file that embeds the compiled SPIR-V shaders. This has to
# happen strictly before the main C++ compiler runs. The binaries that
# use Vulkan have this compile-vk-shaders as a dependency, and split
# actual C++ compiling to a separate subsequent step.
compile-vk-shaders:
	cd vk-shaders && $(MAKE)

myricube-bin: compile-vk-shaders
	$(MAKE) link-myricube-bin

link-myricube-bin:  $(OBJS) $(GL_OBJS) $(VK_OBJS)
	$(CXX)      $(OBJS) $(GL_OBJS) $(VK_OBJS) -ldl -lpthread -lvulkan -o myricube-bin

myricube-vk-bin: compile-vk-shaders
	$(MAKE) link-myricube-vk-bin

link-myricube-vk-bin: $(OBJS) $(GL_DUMMY_OBJS) $(VK_OBJS)
	$(CXX)        $(OBJS) $(GL_DUMMY_OBJS) $(VK_OBJS) -ldl -lpthread -lvulkan -o myricube-vk-bin

myricube-gl-bin: $(OBJS) $(GL_OBJS) $(VK_DUMMY_OBJS)
	$(CXX)   $(OBJS) $(GL_OBJS) $(VK_DUMMY_OBJS) -ldl -lpthread -o myricube-gl-bin

myricube-windows: compile-vk-shaders
	cd windows64 && $(MAKE)



py: libmyricube-cvoxel.so myricube-bin
	echo "NOTE: run './iMyricube.py [your-directory]' to change which world"
	python3 -i iMyricube.py myricube-data/PyWorld/

py-demo: libmyricube-cvoxel.so myricube-bin
	myricube_py_demo=1 python3 iMyricube.py myricube-data/PyDemoWorld/
