///////////////////////////////
     New 3D features of freefem++ 
///////////////////////////////

Version 3.0 is the first version capable of 3D computations.  The price to pay is more difficult mesh generations and difficulties to display the results.  More details can be found in the documentation

1. Mesh generation
==============
There is an interface with the mesh generator tetgen (http://tetgen.berlios.de/)
   see: examples++-3d/Poisson3d.edp  ( mesh of a sphere)
There is a small 3D mesh generator for cylindrical meshes by extrusion of a 2D mesh
   see example: examples++-3d/Lac.edp  for a mesh of a lac
                examples++-3d/Stokes.edp for a mesh a cube.

There is a mesh reader from files in format .mesh of medit (http://www.ann.jussieu.fr/~frey/software.html)
  see example: examples++-3d/dodecaedre01.mesh 
 

2. Graphics
========

2D Graphics
-------------
The standard graphic windows now use a client/server architecture based on openGL GLUT, so the commands have changed:
 - to go to the next graphic type "enter key"
 - to close the graphic window type escape
 - to find out about other commands type the ? key
 
Warning:  you must have "ffglut" in your path; this is done for you by the install program on windows system if you do not uncheck the click box.

3D Graphics
-------------
Real 3D graphics can be obtained with medit, built in freefem andit can be  launched independently of freefem by the command shell "ffmedit".
Freefem can also launch medit by the keywork medit.  Two calls to medit will launch two instances of medit. To close these windows use escape.