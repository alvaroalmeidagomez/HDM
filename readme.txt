===============================================================================
  PROJECT FUNDED BY CENTRO DE MODELAMIENTO MATEMATICO (CMM), SANTIAGO, CHILE
===============================================================================

PROJECT TITLE: 
Numerical Verification of the Heat Equation for 1-Forms on the Unit 2-Sphere

DESCRIPTION:
This experiment consists of verifying the heat equation for 1-forms using the 
Hodge Laplacian operator (Delta_H = d*delta + delta*d) on the two-dimensional 
unit sphere. 

The project contains two main directories, each focusing on the numerical 
solution for a specific initial vector field:

1. Folder V1: Killing Vector Field
   Associated with the numerical solution for the vector field:
   V1(x,y,z) = (-y, x, 0)
   This is a Killing vector field representing a rotation around the z-axis.

2. Folder V2: Projected Vector Field
   Associated with the numerical solution for the vector field:
   V2(x,y,z) = (1,1,1) - <(1,1,1), (x,y,z)> * (x,y,z)
   This represents the orthogonal projection of the constant vector (1,1,1) 
   onto the tangent space at each point of the sphere.

TECHNICAL NOTES:
- Implementation: MATLAB 2017b.
- Architecture: The numerical implementation logic is identical in both 
  folders. The only difference is the vector field definition, which is 
  stored in the "VF.mat" file within each directory.

ACKNOWLEDGEMENTS:
This work was supported by the Centro de Modelamiento Matematico (CMM), 
BASAL fund FB210005 for centers of excellence from ANID-Chile.

===============================================================================
