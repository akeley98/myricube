# Keep It Simple Stupid C/C++/GLSL build dependency tool

cckiss is a tool that automatically tracks the dependencies of C, C++,
and GLSL files and recompiles them if any source or dependency changes
are detected. cckiss is designed to supplement a Unix makefile, rather
than replace `make` entirely. cckiss provides patterns that automate
away the task of writing these kinds of makefile entries:

    bin/foo.o: foo.c foo.h bar.h xyzzy.h
            $(CC) $(CFLAGS) $(CPPFLAGS) foo.c -c -o bin/foo.o

This automation eliminates a whole class of errors related to
forgetting a (possibly indirectly) included header in the list
of dependencies.

Licensed under your choice of CC0 and WTFPL.


# Quick Start

(A more substantial tutorial is provided in `cckiss_EXAMPLE/README.md`)

Clone `cckiss` into the directory where your top-level `Makefile` is.
Include `cckiss/Makefile` in your `Makefile`:

    -include cckiss/Makefile

`cckiss/Makefile` provides patterns that compile source files to
object files named `cckiss/[source-file-path].o` (or `.s`, if you
prefer to compile to assembly), so you only have to provide a rule for
linking all objects together. The variables `CPPFLAGS`, `CC`,
`CFLAGS`, `CXX`, `CXXFLAGS`, and `MAKEFLAGS` are supported in the
standard way.

As an example, if your program contains two source files,
`src/peach.c++` and `src/util.c`, cckiss expects to compile them to
the object files `cckiss/src/peach.c++.o` and `cckiss/src/util.c.o`,
and your `Makefile` should look something like

    default: program

    -include cckiss/Makefile

    CPPFLAGS=-I include
    CC=gcc
    CXX=g++
    CFLAGS=-O2 -Wall
    CXXFLAGS=-O2 -Wall

    # cckiss automatically compiles the listed object/asm files
    program: cckiss/src/peach.c++.o cckiss/src/util.c.o
            $(CXX) cckiss/src/peach.c++.o cckiss/src/util.c.o -o program

The above `Makefile` is the moral equivalent of

    default: program

    bin/peach.o: src/peach.c++ include/peach.h++ some-other-deps.h
            g++ src/peach.c++ -I include -O2 -Wall -c -o bin/peach.o

    bin/util.o: src/util.c include/util.h some-other-deps.h
            gcc src/util.c -I include -O2 -Wall -S -o bin/util.o

    program: bin/peach.o bin/util.o
            g++ bin/peach.o bin/util.o -o program

A few things to note:

1. cckiss supports compiling both C and C++. Any source file ending in
lowercase `.c` is assumed to be C; `.glsl` is glsl, and all others are
assumed C++ (there's no end to the list of C++ file extensions
invented, so I just support them all).

2. cckiss appends `.o` (or `.s`) to the compiled file name, instead of
replacing the extension. So `foo.c` becomes `cckiss/foo.c.o`, NOT
`cckiss/foo.o`.

3. cckiss strictly separates the preprocessing step from the compiling
step. `CPPFLAGS` is only passed to the preprocess stage, while both
`CPPFLAGS` and either `CFLAGS` or `CXXFLAGS` are passed to the C or
C++ compiler stage, respectively. (Note, this currently makes
precompiled headers more-or-less unusable, a significant hole in what
I've set up).

4. You may want to provide a default target before including
`cckiss/Makefile` -- otherwise, the default target for the project
will be the `cckiss` executable itself, which is probably not what you
want.


# Requirements

cckiss requires a Unix-ey environment, and that the C/C++ preprocessor
create directives of the form `# [line number] "included file"`
whenever it includes a file (`#line` is also acceptable). cckiss scans
for these directives in order to know the included dependencies of
each source file. As far as I know, both `gcc` and `clang` do this
correctly.


# Experimental GLSL Support

cckiss supports compiling GLSL files to object (or asm) files that embed
a SPIR-V bytecode array. This is done for any target file of the form

    cckiss/[file-name].glsl.[s|o]

which is compiled from the source file

    [file-name].glsl

As with C files, include files are tracked and recompiling is only
done if any include files (or the source file) has its timestamp
updated.

The glsl file is compiled into a C file declaring a `uint32_t` array
of SPIR-V bytecode, along with a `size_t` declaring the code size in
bytes. The array's name is `[file-name].glsl` with non alphanumeric
characters replaced with `_`; the size is that variable name prefixed
with `sizeof_`. For example, if the shader file is stored in

    vk-shaders/raycast.frag.glsl

Then the shader code and size can be accessed with the declarations

    extern uint32_t vk_shaders_raycast_frag_glsl[];
    extern size_t sizeof_vk_shaders_raycast_frag_glsl;

GLSL compiling is affected by these Makefile variables

`GLSLC`: glsl compiler, for now, assumed to be some version of
`glslangValidator` (`glslc` uses an incompatible CLI).

`GLSLARGS`: arguments passed to `GLSLC`.

`CC`, `CCFLAGS`, `CPPFLAGS`: used for C compiler.


# Motivation

I'm kind of a stubborn Unix hacker and all my C and C++ projects so
far have been built by makefile. I kept thinking all this time that
I've got to find a way to automate the process of listing dependencies
(especially when I make a mistake). `Learn CMake` has been on my todo
list for months now, but so far it just seems far too complicated for
me to tolerate, at least for a personal project without more
experienced people to help me get started with CMake. So, I had to
create an alternative, and here it is.

I'm sure there's some good reason all these industry folks rely on
tools like CMake, but for now, I just don't get it yet. I don't have
to support platforms like Windows and mobile, and I just need
something to help me compile my hobby projects (and who knows, maybe
cckiss is practical enough to be used in production one day, at least
for Unix-only projects).

As for the task itself of scanning for include dependencies: there are
a bunch of tools out there for tracing the inclusion path of a C file,
but that all seemed rather silly to me, since the preprocessor literally
already does that work for us. There's already a ton of directives like

    # 561 "include/mediocre.h"

and we can just piggyback off them to know EXACTLY what headers a
certain source file depends on (even with bizarre nested conditional
includes controlled by preprocessor directives). So, this cckiss tool
(which I've been dreaming of for years) turned out to be really simple
to write: all I had to do was preprocess the file, grep out a list of
header dependencies, and check for modifications in the listed files
to know if a recompilation is needed.  Stupid Simple, what was I
waiting for?


# Implementation

Everything is powered by the `cckiss/cckiss` executable, built by the
provided makefile from `cckiss/cckiss.cc` (C is for cookie...) This
executable is targetted by C/C++ patterns like:

    cckiss/%.c.s : .cckiss.PHONY cckiss/cckiss
	    cckiss/cckiss $@ $(CC) .cckiss.CPPFLAGS $(CPPFLAGS) .cckiss.CXXFLAGS $(CFLAGS)

When called for the first time on a target file, `cckiss/cckiss` strips away
the `cckiss/` prefix and `.o` (or `.s`) suffix to recover the path to
the source file to be compiled. It then

1. Preprocesses the source file, storing it in
`cckiss/[path-to-source].i` (or `.ii`, for C++ files). This step uses
the makefile variables `CPPFLAGS`, `CC`, and `CFLAGS` for C;
`CPPFLAGS`, `CXX`, and `CXXFLAGS` for C++ files.

2. Scans the preprocessed file for a list of dependency files (using the
`# [line-number] [file-name]` directives). This list of dependency
files is stored in `cckiss/[path-to-source]-deps.txt`. Each file listed
is separated with exactly 1 newline, with exactly 0 or 1 trailing newlines
in the file. I considered it unlikely that any source file path would contain
a newline (how could such a file be included under C syntax rules?) However,
the newlines could be replaced with `'\0'` if this is a real problem.

3. Compile the preprocessed source file to create the target object or
assembly file. This uses the makefile variables `CC` and `CFLAGS` for
C files, `CXX` and `CXXFLAGS` for C++.

4. Future invokations of `cckiss/cckiss` on the same target file look
up the `-deps.txt` file, and skip recompilation if all listed files
have modification times before the modification time of the target file.
However, if `B` is in the environment variable `MAKEFLAGS`, recompilation
is done unconditionally, for consistency with `make`.


# TODO (consider)

The `-deps.txt` files include a lot of system headers. Consider providing an
option that allows skipping checks for such dependencies.

There's no way to tailor `CFLAGS` and `CXXFLAGS` per-file. This could
be useful, for example, to provide flags like `-Ofast` or
`-funroll-loops` that may be needed for certain performance-critical
files but are otherwise unsafe for general use. Could this be done by
allowing meta-data in source files (e.g. a `CCKISS_ARG` preprocessor
macro), or by providing a separate args file (e.g. `foo.c.cckiss` for
`foo.c`)?

The `-MD` option for `gcc` and `clang` can generate a list of
dependency files, making the separate preprocess stage
unneccessary. However, this is still needed for `glslangValidator`,
which does not (as far as I can tell) support such a feature. `glslc`
claims to support this flag, but it does not support compiling to a C
file (and also does not recognize compound file extensions,
e.g. `.frag.glsl`, a gratuitious incompatibility with
`glslangValidator`).
