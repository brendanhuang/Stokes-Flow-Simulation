# Stokes-Flow-Simulation
Stokes-Flow-Simulation is a Matlab implementation of boundary element method (BEM) and method of fundamental solutions (MFS) for simulation of Stokes flow based on traction and velocity boundary conditions.

<img src="/images/flowfield_streamfunction.png" width = "40%" align="middle">

This repository contains an implementation of numerical simulations of low Reynold's number flow (Stokes flow). The work was done as part of my doctoral thesis at Yale University [1]. There are three possible types of simulations that can be performed with the code:

<ol>
<li> The method of fundamental solutions (MFS) to solve for flow in 2D
<li> The boundary element method (BEM) to solve for flow in 2D
<li> The BEM to solve for flow in 3D
</ol>

(MFS in 3D has been left out due to the unregularized version of MFS being an ill-posed problem and thus sensitive to numerical errors.) In all cases, the numerics solve for the flow in the domain after specifying traction and flow boundary conditions. The default setting is to simulate a geometry similar to a [lid-driven cavity](https://www.cfd-online.com/Wiki/Lid-driven_cavity_problem).

# Installation

1. Download folder containing m-files.
2. Add all folders and subfolders to path in Matlab.
3. Open doit_sim_BEM_2D.m and execute cell by cell.

# How to use this repository
The repository contains a series of .m files as well as a tutorial document (StokesFlowSimulation_Tutorial.pdf). The m-files are separated into "doit" files that can be run immediately. These doit files call back end functions.  All of the functions are explained in the tutorial document. 

The tutorial explains the underlying theory, technical details, and structure of the implementation of these numerical simulations. The tutorial itself is an appendix in my thesis, the full version of which can also be downloaded on my github page (@brendanhuang).

# References
A fair amount of the code here is adapted from [2], an implementation of BEM for Laplace's equation. The Stokes flow implementation is follows the overall scheme of [3]. Additional background can be found in [] and [].

[1] B.K. Huang. All optical quantification of ciliary physiology. PhD thesis, Yale University, 2015.

[2] Stephen Kirkup, Javad Yazdani, NE Mastorakis, M Poulos, V Mladenov, Z Bo- jkovic, D Simian, S Kartalopoulos, A Varonides, and C Udriste. A gentle introduction to the boundary element method in matlab/freemat. In WSEAS International Conference. Proceedings. Mathematics and Computers in Science and Engineering, number 10. WSEAS, 2008.

[3] C. Pozrikidis. A practical guide to boundary element methods with the software library BEMLIB. Chapman & Hall/CRC, Boca Raton, 2002.

[4] C. Pozrikidis. Boundary integral and singularity methods for linearized viscous flow. Cambridge University Press, Cambridge England; New York, 1992.

[5] C. Pozrikidis. Introduction to theoretical and computational fluid dynamics. Oxford University Press, Oxford England; New York, 1997.
