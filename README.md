# Stokes-Flow-Simulation
A Matlab implementation of boundary element method (BEM) and method of fundamental solutions (MFS) for simulation of Stokes flow based on traction and velocity boundary conditions.

This repository contains an implementation of numerical simulations of low Reynold's number flow (Stokes flow). The work was done as part of my doctoral thesis at Yale University. There are three possible types of simulations that can be performed with the code 

<ol>
<li> The method of fundamental solutions (MFS) to solve for flow in 2D
<li> The boundary element method (BEM) to solve for flow in 2D
<li> The BEM to solve for flow in 3D
</ol>

(MFS in 3D has been left out due to the unregularized version of MFS being a Hadamard ill-posed problem.) In all cases, the numerics solve for the flow in free space after specifying traction and flow boundary conditions.

# How to use this repository
The repository contains a series of .m files.  Backend files . There are several examples files ().

There is also a tutorial explaining the technical details on the implementation of these numerical simulations. The tutorial is an appendix in my thesis, the full version of which can also be downloaded on this website.

# References
A fair amount of the code here is adapted from [], an implementation of BEM for Poisson's equation. The Stokes flow implementation is follows the overall scheme of []. Additional background can be found in [] and [].

# Tutorial files:

### Testing

<img src="/images/maxresdefault.jpg" width = "40%">
