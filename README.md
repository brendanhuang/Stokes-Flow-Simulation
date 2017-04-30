# Stokes-Flow-Simulation
Stokes-Flow-Simulation is a Matlab implementation of boundary element method (BEM) and method of fundamental solutions (MFS) for simulation of Stokes flow based on traction and velocity boundary conditions.

<img src="/images/flowfield_streamfunction.png" width = "80%" align="middle">

This repository contains an implementation of numerical simulations of low Reynold's number flow (Stokes flow). The work was done as part of my doctoral thesis at Yale University [1]. There are three possible types of simulations that can be performed with the code:

<ol>
<li> The method of fundamental solutions (MFS) to solve for flow in 2D
<li> The boundary element method (BEM) to solve for flow in 2D
<li> The BEM to solve for flow in 3D
</ol>

In all cases, the routines numerically solve for the vectorial flow field in the interior of the domain after specifying traction and/or flow boundary conditions. The default setting is to simulate a geometry similar to a [lid-driven cavity](https://www.cfd-online.com/Wiki/Lid-driven_cavity_problem). In certain cases, the pressure field, shear stress tensor, and/or stream function can also be directly computed.

# Installation

1. Download folder containing m-files.
2. Add all folders and subfolders to path in Matlab.
3. Open doit_sim_BEM_2D.m and execute cell by cell.

# How to use this repository
The repository contains a series of m-files as well as a tutorial document. The m-files in turn are separated into "doit" executable files that can be run immediately. These files are all located in the scripts folder. The executable files in turn call back-end functions. The functions are divided folderwise into bem_2d_functions, bem_3d_functions, and mfs_2d_functions based on which simulation they are called by.  In addition, there is a common shared folder of numerical routines used by all simulation methods, as well as some visualization scripts. 

The tutorial (StokesFlowSimulation_Tutorial.pdf) explains the underlying theory, technical details, and structure of the implementation of these numerical simulations. All of the scripts included are explained in the tutorial document. The tutorial itself is an appendix in my thesis, the full version of which can also be downloaded on my github page (@brendanhuang).

# References
A couple subroutines are adapted from [2], an implementation of BEM for Laplace's equation. The Stokes flow implementation follows the overall scheme of [3]. Additional background can be found in [4] and [5].

[1] B.K. Huang. All optical quantification of ciliary physiology. PhD thesis, Yale University, 2015.

[2] Stephen Kirkup, Javad Yazdani, NE Mastorakis, M Poulos, V Mladenov, Z Bo- jkovic, D Simian, S Kartalopoulos, A Varonides, and C Udriste. A gentle introduction to the boundary element method in matlab/freemat. In WSEAS International Conference. Proceedings. Mathematics and Computers in Science and Engineering, number 10. WSEAS, 2008.

[3] C. Pozrikidis. A practical guide to boundary element methods with the software library BEMLIB. Chapman & Hall/CRC, Boca Raton, 2002.

[4] C. Pozrikidis. Boundary integral and singularity methods for linearized viscous flow. Cambridge University Press, Cambridge England; New York, 1992.

[5] C. Pozrikidis. Introduction to theoretical and computational fluid dynamics. Oxford University Press, Oxford England; New York, 1997.
