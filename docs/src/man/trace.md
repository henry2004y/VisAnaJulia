# Streamline and Trajectory

Streamline and trajectory are related topics in physical modeling, common seen in fluid and particle simulations.
In the future, the stream tracing should become a stand-alone package!

## Streamline Tracing

Given an unstructured grid with node points and connectivity, how should you do the streamline tracing?

Brute force algorithm:
1. find the grid cell you are currently in;
2. move along the vector direction until you hit the boundary of that cell;
3. find the neighbour who shares the same edge;
4. use the vector direction in the next cell and move along that direction;
5. repeat 2-4 until you hit the preset stopping criteria.

Some questions during the process:
* How to find the neighbouring cell?
* How to determine which boundary edge will you cross?
* How to improve the search speed?
* How to improve accuracy?

First make it work, then make it better and fast.

I found an approach called Pollock method.

I need an adaptive step control integration scheme like rk45.

There is an implementation of streamline tracing in Matlab called tristream.

There is another implementation in yt library:
Streamlining through a volume is useful for a variety of analysis tasks. By specifying a set of starting positions, the user is returned a set of 3D positions that can, in turn, be used to visualize the 3D path of the streamlines. Additionally, individual streamlines can be converted into `YTStreamline` objects, and queried for all the available fields along the streamline.

The implementation of streamlining in yt is described below.

1. Decompose the volume into a set of non-overlapping, fully domain tiling bricks, using the `AMRKDTree` homogenized volume.

2. For every streamline starting position:

  * While the length of the streamline is less than the requested length:

    1. Find the brick that contains the current position

    2. If not already present, generate vertex-centered data for the vector fields defining the streamline.

    3. While inside the brick

      1. Integrate the streamline path using a Runge-Kutta 4th order method and the vertex centered data.

      2. During the intermediate steps of each RK4 step, if the position is updated to outside the current brick, interrupt the integration and locate a new brick at the intermediate position.

3. The set of streamline positions are stored in the `Streamlines` object.

## Particle Tracing

I have a plan of incorporating particle tracing into this module. WIP
