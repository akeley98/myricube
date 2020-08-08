default: run

CPPFLAGS=-I glfw/include -I ./glad/include -I .
CXXFLAGS=-Wall -Wextra -O3 -g -std=c++17
CFLAGS=  -Wall -Wextra -O3 -g

-include cckiss/Makefile

glfw-build/Makefile:
	cd glfw-build && cmake ../glfw

glfw-build/src/libglfw3.a: glfw-build/Makefile
	cd glfw-build && make

run: myricube-planet-bin
	./myricube-planet

all: myricube-planet-bin \
     myricube-axis-bin \
     cckiss/apps/axis.cc.s \
     myricube-congestion-bin \
     cckiss/apps/congestion.cc.s \
     myricube-random-walk-bin \
     cckiss/apps/random-walk.cc.s \
     myricube-skygrid-bin \
     cckiss/apps/skygrid.cc.s \
     myricube-hexload-app-bin \
     cckiss/apps/hexload-app.cc.s \
# Empty line for backslash

OBJS=cckiss/myricube.cc.s \
     cckiss/window.cc.s \
     cckiss/renderer.cc.s \
     cckiss/glad/src/glad.c.s \
     glfw-build/src/libglfw3.a \
# Empty line for backslash

myricube-planet-bin: $(OBJS) cckiss/apps/marlo-planet.cc.s
	$(CXX) $(OBJS) cckiss/apps/marlo-planet.cc.s -ldl -lGL -lpthread -o myricube-planet-bin


myricube-%-bin: $(OBJS) cckiss/apps/%.cc.s
	$(CXX) $(OBJS) cckiss/apps/$*.cc.s -ldl -lGL -lpthread -o myricube-$*-bin
