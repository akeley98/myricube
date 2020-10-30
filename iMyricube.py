#! /usr/bin/python3 -i
# Hacked-together Python client for myricube.

import sys
import os
import ctypes
file_dir = os.path.split(__file__)[0]
lib = ctypes.cdll.LoadLibrary(os.path.join(file_dir, "./libmyricube-cvoxel.so"))

if len(sys.argv) != 2: raise Exception("Need 1 cmd line arg: world filename")
world_filename = sys.argv[1]

code = lib.myricube_select_world(bytes(world_filename, "utf-8"))
if code != 0: raise Exception("myricube_select_world failed")
sys.stderr.write("Opened world '%s'\n" % world_filename)

"""Set the voxel at the (x,y,z) coordinate of the current world to the
given color (8-bit red, green, blue)"""
def set_rgb(x, y, z, r, g, b):
    lib.myricube_set(int(x), int(y), int(z), rgb(r, g, b))

"""Convert 8-bit RGB triple to 32-bit opaque voxel value (last argument to Set)"""
def rgb(r, g, b):
    return min(max(r, 0), 255) << 24 | min(max(g, 0), 255) << 16 | min(max(b, 0), 255) << 8 | 128

"""Set the voxel at the (x,y,z) coordinate of the current world to the
given 32-bit voxel value. NOTE: 0 is an empty voxel."""
def Set(x, y, z, voxel32):
    lib.myricube_set(int(x), int(y), int(z), int(voxel32))

"""Get the 32-bit value of the voxel at the given coordinate of the
current world. NOTE: At the moment, the underlying C library may
create a new chunk group on disk, so this is not a truly read-only
operation."""
def get(x, y, z):
    return lib.myricube_get(int(x), int(y), int(z))

"""Fill the box with corners (x0, y0, z0) and (x1, y1, z1) (inclusive)
with the given voxel value."""
def fill(x0, y0, z0, x1, y1, z1, voxel32):
    return lib.myricube_fill(
        int(x0), int(y0), int(z0), int(x1), int(y1), int(z1), int(voxel32))

def zholes(x0, y0, z0, x1, y1, z1, voxel32_0, voxel32_1):
    return lib.myricube_zholes(
        int(x0), int(y0), int(z0), int(x1), int(y1), int(z1), int(voxel32_0), int(voxel32_1))

# Launch myricube to view the current world.
os.environ["myricube_app"] = "ViewWorld"
os.environ["myricube_world"] = world_filename
pid = os.fork()
if pid == 0:
    myricube_bin = os.path.join(file_dir, "./myricube-bin")
    os.execvp(myricube_bin, (myricube_bin,))

# Launch demo if myricube_py_demo environment variable is set to nonzero value.
do_demo = False
try:
    do_demo = int(os.environ["myricube_py_demo"])
except KeyError:
    pass

# Basically the pipes screensaver.
if do_demo:
    from random import randrange
    from random import choice
    from time import sleep
    radius = 133

    up = 0
    down = 1
    left = 2
    right = 3
    forward = 4
    back = 5

    while 1:
        # Clear the world periodically.
        fill(-radius, -radius, -radius, radius, radius, radius, 0)

        # One pipe per iteration
        for pipes in range(18):
            low = -radius // 2
            high = +radius // 2
            x = randrange(low, high)
            y = randrange(low, high)
            z = randrange(low, high)

            # Pick a color for this pipe.
            n = randrange(3)

            if n == 0:
                r = randrange(88, 255)
                g = randrange(88, 255)
                b = 88
            elif n == 1:
                r = randrange(88, 255)
                g = 88
                b = randrange(88, 255)
            else:
                r = 88
                g = randrange(88, 255)
                b = randrange(88, 255)

            direction = right
            pipe_in_bounds = True

            # One straight segment of the pipe per iteration.
            while pipe_in_bounds:
                choices = (
                    (left, right, forward, back),
                    (left, right, forward, back),
                    (up, down, forward, back),
                    (up, down, forward, back),
                    (up, down, left, right),
                    (up, down, left, right))[direction]

                direction = choice(choices)

                # Generate one straight segment of pipe in the chosen direction
                for iter in range(randrange(10, 20)):
                    # The money function
                    set_rgb(x, y, z, r, g, b)

                    if direction == up: y += 1
                    elif direction == down: y -= 1
                    elif direction == left: x -= 1
                    elif direction == right: x += 1
                    elif direction == forward: z += 1
                    else: z -= 1

                    pipe_in_bounds &= x >= -radius and x <= radius
                    pipe_in_bounds &= y >= -radius and y <= radius
                    pipe_in_bounds &= z >= -radius and z <= radius
                    sleep(0.001)
                    if not pipe_in_bounds: break

                    # You don't have to worry about this.  This is
                    # just to get the demo to exit cleanly when the
                    # myricube window is closed.
                    pid, status = os.waitpid(pid, os.WNOHANG)
                    if pid != 0: sys.exit(0)
