# Makefile for converting glsl shaders to SPIR-V, then to .xxd files
# (lists of comma-separated byte values).

# Compile unconditionally using phony target due to no way to track
# include dependencies.
%.spv.xxd: .vk-shaders.PHONY
	glslangValidator $* -V -o $*.spv
	xxd -i < $*.spv > $@
# glslangValidator {glsl input file} {-V:actually compile} -o {spv file}
# xxd {-i:as include file} < {spv file} > {.xxd include file}

.vk-shaders.PHONY :
	@echo -n ""

default: background.vert.spv.xxd background.frag.spv.xxd mesh.vert.spv.xxd mesh.frag.spv.xxd raycast.vert.spv.xxd raycast.frag.spv.xxd
