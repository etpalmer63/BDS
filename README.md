# BDS


Start by using HeatEqn 0 as template. 

Then inital problem values are in the 3D paper.

replace the evolve section of the code with call to bds.


assume periodic boundary conditions 

s should need 3-4 ghost cells.

Will need to add multifab to hold all the u,v,w MAC velocities.

For the fortran `allocate` command, will need to pass in a tmp multifab. 
(Can't be allocated on the fly)


A big part of porting the code will be to figure out which kind
of box is needed for each set of nested for loops. For example,
-1 to +1 is cell-centered with 1 ghost cell
-1 to +2 is nodal with 1 ghost cell


s_old will likely need 4 ghost cells
