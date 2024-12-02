---
name: crack-3d
category: Solid Mechanics
folder: 3d
---

## Build a 3d mesh with a fissure, check eps in trunc check mesh
We begin with a crack in a circle
~~~freefem
real theta = 1*pi/180.; // 0.7 degree  limit to build the 2d mesh if no refinement in zero
real xt = cos(theta),yt=sin(theta);
cout << xt << " " << yt << endl;
border C(t=theta, 2*pi-theta) { x = cos(t); y=sin(t); label=1;}
border Ss(t=0,1;i){ t= t; x= t*xt; y = t*yt;if(i) y = -y; label=2+i;}
real hs = 0.05;
int nC = (2*pi-2*theta)/hs;
int ns = (1./hs)/1;//   
int[int] nS=[ ns,-ns ];
func bord = C(nC)+Ss(nS); 
plot(bord,wait=1); // check the border is closed
mesh Th2 = buildmesh(bord,nbvx=1000000);
plot(Th2);
~~~

![][_crack2d]

Let us verify that the mesh points are not closed to one another
~~~freefem
int nt = Th2.nt, nv = Th2.nv;
Th2=trunc(Th2,true);// verif to closeness of points
assert(nt==Th2.nt);
assert(nv==Th2.nv);
~~~
The 3d mesh is obtained by a vertical extrusion using $\texttt{buildlayers}$
~~~freefem
mesh3 Th = buildlayers(Th2,1./hs);
plot(Th,wait=1);
 nt = Th.nt;
 nv = Th.nv;
// check mesh
Th=checkmesh(Th); // works only in 3d hence not used above in 2d
assert(nt==Th.nt);
assert(nv==Th.nv);
cout<<"All is well that ends well "<<endl;
~~~

| The final result: |
|-------------------|
|![][_crack3d]      |


[_crack2d]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/crack-3d/crack2d.png

[_crack3d]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/crack-3d/crack3d.png
