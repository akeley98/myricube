default: all myricube-random-walk-bin
	./myricube-random-walk

CPPFLAGS=-I glfw/include -I ./glad/include -I include -I glm
CXXFLAGS=-Wall -Wextra -O3 -g -std=c++17
CFLAGS=  -Wall -Wextra -O3 -g

-include cckiss/Makefile
-include OBJS-list

glfw-cmake:
	cd glfw-build && cmake ../glfw

glfw-build/src/libglfw3.a: glfw-cmake
	cd glfw-build && make

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

# Target all-w64 to REALLY make everything for Windows and Linux.
all-w64: all
	cd windows64 && make

myricube-planet-bin: $(OBJS) cckiss/apps/marlo-planet.cc.s
	$(CXX) $(OBJS) cckiss/apps/marlo-planet.cc.s -ldl -lGL -lpthread -o myricube-planet-bin


myricube-%-bin: $(OBJS) cckiss/apps/%.cc.s
	$(CXX) $(OBJS) cckiss/apps/$*.cc.s -ldl -lGL -lpthread -o myricube-$*-bin
