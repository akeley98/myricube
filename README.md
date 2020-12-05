![myricube screenshot](./screenshot.png)

## Myricube – Hybrid Mesh/Raycasting Voxel Renderer

Welcome to the rudimentary readme. This project is a work-in-progress.
A technical essay on my goals in this project is in `src/myricube.cc`.

## Windows 64-bit

You need OpenGL 4.5+ or Vulkan 1.1+ graphics drivers for this to
work. Download the project as a zip, unzip it, and run the precompiled
exe file, either

`myricube-master/myricube-master/windows64/myricube.exe` for Vulkan 1.1
(recommended) or

`myricube-master/myricube-master/windows64/myricube.exe` for OpenGL 4.5

I have basically zero experience with Windows programming, so if the
exe doesn't work and you feel nice, file an issue giving a description
of what error you see. At the moment, I don't support building on
Windows. I use several Unix-native tools like Makefiles and `xxd`
(to convert SPIR-V shaders to C header files) that I don't feel
like porting to Windows. Everything is cross-compiled on Linux for Windows.

## GNU / Linux

You need OpenGL 4.5 and `git`. For my modified version of
GLFW, you also need `cmake` and `xorg-dev` (not sure what you need
if you use Wayland).

Git clone the repo and run `make -j` to build and run the OpenGL
version. If you don't have the GLM library you have to clone this
repository recursively.

If you want to try out the Vulkan version, you need to have the
Vulkan SDK installed in order to compile shaders to SPIR-V. Then,
make `myricube-bin` or `myricube-vk-bin`. This is a bit unfair
as the Windows version just has the shaders baked into the executable,
but, since I'm doing actual development on Linux, I don't want
intermediate files polluting the repository so too bad.

## Controls

If you don't like these controls, look at `myricube-data/default-keybinds.txt`.
In any event, don't get too used to these as I'll probably mess everything up
if I ever add a real UI.

`wasd` to move; `q` to slow down, `e` to speed up.

right-click-drag or trackpad two-finger scroll to look around.

left/right arrows or thumb buttons to go through camera history.

`k` in the default app generates another random walk.

`f9`/`f10` to change the render distance (be careful not to run out of memory),
and `f`,`b` to toggle fog/fog color.

## World Files

Myricube supports loading and storing voxel models (worlds) from disk.
To do this, run `./myricube` (Linux) or one of the myricube `.exe` files
with the name of a directory as the argument. Make sure the directory
name ends with a trailing `/` (or `\\` on Windows). This creates a
`world.myricube` file in the directory (if it already exists, it must
either be empty or have been created earlier by Myricube).

On Windows, you can also associate `.myricube` files with `myricube.exe`.
Create new worlds by creating a new folder with an empty `world.myricube`.

This is a little underwhelming for now as I haven't added a UI for
Myricube. For now, I have a Python API (up next); you can use
`./iMyricube.py` with a directory argument in order to open a world
with a Python command line to edit them and a Myricube window to view
them.

## Python API

Myricube uses non-exclusive file memory-mapping to load worlds; thus,
you can modify worlds with external applications and see changes
in real time. I have provided a Python script for doing so: the
aforementioned `iMyricube.py`. You can also use this to create
a voxel animation with Python. See everything under the `if do_demo:`
line of `iMyricube.py` to see an example of this.

Run this demo with `make py-demo` (Linux) or `windows64/pipes.bat`
(Windows). This demo is inspired by the old Windows pipes screensaver.

## C++ Apps

The older way to play with Myricube is to write C++ apps, which are
bits of code that generate an alien planet or simulate cars looping on
the surface of a cube forever or otherwise mess around with the voxel
world. This doesn't work so good on Windows.

Look at `include/app.hh` and `OBJS-list` to start. You can select an
app using the `myricube_app` environment variable,
e.g. `myricube_app=MarloPlanet ./myricube` generates my <!--
brilliant, knockout GORGEOUS --> roommate Marlon's awesome Perlin
noise cave planet. Look in the `apps` directory for more choices.

## Environment Variables

For now you can change the behavior of Myricube with environment
variables. This is kinda hacky and subject to change and out-of-date
documentation.

`myricube_api` selects the graphics API, either `gl` for OpenGL 4.5
core profile or `vk` for Vulkan 1.1. The Vulkan version seems to run
faster but it's not 100% production quality, e.g. it falls flat on its
face if an image format is missing, or if it runs out of GPU
memory. Try the OpenGL version if the Vulkan version doesn't work.

`myricube_app` selects which C++ app to run. By default, it runs the
random walk app if no command line argument to Myricube is provided,
or the `ViewWorld` app if one command line argument is provided.
Interesting choices include `MarloPlanet`, Marlon's cave planet, which
is generated in a background thread so you can see it built in real
time it's SO COOL! Also `Congestion` is neat, it's a cellular automata
game run on the surface of a 3D cube based on the
Biham–Middleton–Levine traffic model.

`myricube_validation` determines whether Vulkan Validation Layers or
OpenGL `glGetError` is enabled. `0` disables (default), nonzero
enables.

`myricube_endian_fix` is a real hack. The file format for worlds changed
endian-ness through development, setting this environment variable to
a nonzero value enables auto-fixing the old format to new.

`myricube_world` determines the file to load for the `ViewWorld` app.
You probably don't have to worry about this.
