cckiss/cckiss : cckiss/cckiss.cc
	$(CXX) -std=c++14 -O2 cckiss/cckiss.cc -o cckiss/cckiss -Wall -Wextra -Wno-constant-logical-operand -Wno-unused-function

.cckiss.PHONY :
	@echo -n ""

# Provide default glsl compiler.  For now, only glslangValidator is
# supported (but you may still need to override to a specific path or
# location)
ifeq ($(origin GLSLC),undefined)
GLSLC := glslangValidator
endif

# Compile .glsl to SPIR-V, then to C object file.
cckiss/%.glsl.s : .cckiss.PHONY cckiss/cckiss
	cckiss/cckiss $@ $(CC) .cckiss.CPPFLAGS $(CPPFLAGS) .cckiss.CXXFLAGS $(CFLAGS) .cckiss.GLSLC $(GLSLC) .cckiss.GLSLARGS $(GLSLARGS)

cckiss/%.glsl.o : .cckiss.PHONY cckiss/cckiss
	cckiss/cckiss $@ $(CC) .cckiss.CPPFLAGS $(CPPFLAGS) .cckiss.CXXFLAGS $(CFLAGS) .cckiss.GLSLC $(GLSLC) .cckiss.GLSLARGS $(GLSLARGS)

# Compile .c to C object file.
cckiss/%.c.s : .cckiss.PHONY cckiss/cckiss
	cckiss/cckiss $@ $(CC) .cckiss.CPPFLAGS $(CPPFLAGS) .cckiss.CXXFLAGS $(CFLAGS)

cckiss/%.c.o : .cckiss.PHONY cckiss/cckiss
	cckiss/cckiss $@ $(CC) .cckiss.CPPFLAGS $(CPPFLAGS) .cckiss.CXXFLAGS $(CFLAGS)

# Other extensions are assumed to be C++, compile to C++ object file.
cckiss/%.s : .cckiss.PHONY cckiss/cckiss
	cckiss/cckiss $@ $(CXX) .cckiss.CPPFLAGS $(CPPFLAGS) .cckiss.CXXFLAGS $(CXXFLAGS)

cckiss/%.o : .cckiss.PHONY cckiss/cckiss
	cckiss/cckiss $@ $(CXX) .cckiss.CPPFLAGS $(CPPFLAGS) .cckiss.CXXFLAGS $(CXXFLAGS)

# NOTE: If you mess with the above file extension matches, you have to
# modify the cckiss executable to match as well.

# Annoying thing needed to avoid Makefile making itself.
cckiss/Makefile:
	@echo -n ""
