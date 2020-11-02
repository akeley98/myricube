![myricube screenshot](./screenshot.png)

## Myricube – OpenGL Hybrid Mesh/Raycasting Voxel Renderer

Welcome to the rudimentary readme. This project is a work-in-progress.
A technical essay on my goals in this project is in `src/myricube.cc`.

## Windows 64-bit

Not supported on Windows at the moment.

<!--You need OpenGL 4.5 graphics drivers for this to work. Download the
project as a zip, unzip it, and run the precompiled exe at
`myricube-master/myricube-master/windows64/myricube.exe`.-->

I have basically zero experience with Windows programming, so if the
exe doesn't work and you feel nice, file an issue giving a description
of what error you see.

## GNU / Linux

You need OpenGL 4.5 and `git`. For my modified version of
GLFW, you also need `cmake` and `xorg-dev` (not sure about Wayland).

Git clone the repo and run `make -j`. If you don't have the GLM
library you have to clone recursively.

## Controls

If you don't like these controls, look at `myricube-data/default-keybinds.txt`.
In any event, don't get too used to these as I'll probably mess everything up
if I ever add a real UI.

`wasd` to move.

right-click-drag or trackpad two-finger scroll to look around.

left/right arrows or thumb buttons to go through camera history.

`k` in the default app generates another random walk.

`f9`/`f10` to change the render distance (be careful not to run out of memory),
and `f`,`b` to toggle fog/fog color.

`f3`/`f4` to change the maximum resolution (not implemented in current
branch). NOTE: Pressing `f4` a few times can work-around weird
z-fighting artifacts (sometimes happens on Windows).

## Playing with Myricube

One day I may add a real UI or maybe a Python API, but for now you have
to be a programmer and create C++ apps, which are bits of code that
generate an alien planet or simulate cars looping on the surface of
a cube forever or otherwise mess around with the voxel world.

Look at `include/app.hh` and `OBJS-list`. You can select an app using
the `myricube_app` environment variable,
e.g. `myricube_app=MarloPlanet ./myricube` generates my <!--
brilliant, knockout GORGEOUS --> roommate Marlon's awesome Perlin
noise cave planet. Look in the `apps` directory for more choices.
