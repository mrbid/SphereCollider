# SphereCollider
Spheres that bounce around.

Originally made for https://github.com/mrbid/ChaoticSystemModelling but seperated out and improved for a bit of fun.

## Controls
- O = Orbit in.
- N = New simulation.
- F = FPS to console.

## Compile
```
sudo apt install libglfw3 libglfw3-dev
gcc main.c glad_gl.c -I inc -Ofast -lglfw -lm -o SphereCollider
```
