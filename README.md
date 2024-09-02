# Summer_Project_2024
Investigation of Regularization by Discretization for the Reconstruction of Randomly Rough Surfaces

This repository contains all of the functions/scripts that were used to numerically implement and investigate the direct and inverse wave scattering problems as part of a summer research project under Dr. Orsola Rath Spivack. 

In the direct problem, a Gaussian beam impinges upon a 1-dimensional randomly rough surface and the values of the surface derivative and scattered field at a height z=Z above the mean surface are calculated at nodes along the surface. The inverse problem then uses the scattered field data generated in the direct problem to successively reconstruct the surface derivative and height data of the surface in a marching fashion. Dirichlet boundary conditions are imposed along the surface. 
