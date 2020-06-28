default: run

CPPFLAGS=-I ./glad/include -I .
CXXFLAGS=-Wall -Wextra -O3 -g -std=c++17
CFLAGS=  -Wall -Wextra -O3 -g

-include cckiss/Makefile

run: myricube-planet-bin
	./myricube-planet

all: myricube-planet-bin \
     myricube-axis-bin \
     myricube-congestion-bin \
     myricube-random-walk-bin \
     myricube-skygrid-bin \
     myricube-hexload-app-bin \
# Empty line for backslash

OBJS=cckiss/myricube.cc.s \
     cckiss/window.cc.s \
     cckiss/renderer.cc.s \
     cckiss/glad/src/glad.c.s \
# Empty line for backslash

myricube-planet-bin: $(OBJS) cckiss/apps/marlo-planet.cc.s
	$(CXX) $(OBJS) cckiss/apps/marlo-planet.cc.s -ldl -lGL -lSDL2 -o myricube-planet-bin


myricube-%-bin: $(OBJS) cckiss/apps/%.cc.s
	$(CXX) $(OBJS) cckiss/apps/$*.cc.s -ldl -lGL -lSDL2 -o myricube-$*-bin
