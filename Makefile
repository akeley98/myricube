default: myricube-gl-bin
	./myricube

CPPFLAGS=-I glfw/include -I ./glad/include -I include -I glm -I vk-helpers -I vk-shaders
CFLAGS=  -Wall -Wextra -Wno-missing-field-initializers -O3 -g -fPIC
CXXFLAGS=$(CFLAGS) -std=c++17

-include cckiss/Makefile
-include OBJS-list

glfw-build/Makefile:
	cd glfw-build && cmake ../glfw >../build/glfw-cmake.stdout

do-glfw-build: glfw-build/Makefile
	cd glfw-build && $(MAKE) >../build/glfw-make.stdout 2>../build/glfw-make.stderr

libmyricube-cvoxel.so: $(LIBOBJS)
	$(CXX) $(LIBOBJS) -shared -o libmyricube-cvoxel.so

all: myricube-bin myricube-gl-bin myricube-vk-bin myricube-windows libmyricube-cvoxel.so

myricube-bin:  $(OBJS) $(GL_OBJS) $(VK_OBJS)
	$(CXX) $(OBJS) $(GL_OBJS) $(VK_OBJS) -ldl -lpthread -lvulkan -o myricube-bin

myricube-vk-bin: $(OBJS) $(GL_DUMMY_OBJS) $(VK_OBJS)
	$(CXX)   $(OBJS) $(GL_DUMMY_OBJS) $(VK_OBJS) -ldl -lpthread -lvulkan -o myricube-vk-bin

myricube-gl-bin: $(OBJS) $(GL_OBJS) $(VK_DUMMY_OBJS)
	$(CXX)   $(OBJS) $(GL_OBJS) $(VK_DUMMY_OBJS) -ldl -lpthread -o myricube-gl-bin

myricube-ninja: do-glfw-build
	build/generate_ninja.py > build.ninja && \
	ninja -v -t compdb cc cxx link > compile_commands.json

myricube-windows:
	cd windows64 && $(MAKE)



py: libmyricube-cvoxel.so myricube-gl-bin
	echo "NOTE: run './iMyricube.py [your-directory]' to change which world"
	python3 -i iMyricube.py myricube-data/PyWorld/

py-demo: libmyricube-cvoxel.so myricube-gl-bin
	myricube_py_demo=1 python3 iMyricube.py myricube-data/PyDemoWorld/
