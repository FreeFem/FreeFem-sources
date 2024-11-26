---
name: extract-boundary-3d
category: mesh
folder: 3d
---

## Extract from the mesh the boundary with a given label

The cube has 6 faces, each with a different label. The function $\texttt{extract}$ returns a surface mesh on which all vertices have the given label
~~~freefem
load "medit"
int n= 10;
int nvb = (n+1)^3 - (n-1)^3;// Nb boundary vertices
int ntb = n*n*12; // Nb of Boundary triangle 
mesh3 Th=cube(n,n,n); // lalels are from 1 to 6
int[int] ll=[1,3];//for example
meshS ThS=extract(Th,label=ll);// extract boundary of 3d Mesh with given label
plot(ThS);
~~~

| The surface mesh extracted |
|----------------------------|
|![][_solution]              |

[_solution]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/extract-boundary3d/solution.png