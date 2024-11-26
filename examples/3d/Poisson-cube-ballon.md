---
name: Periodic-cube-ballon
category: math
layout: 3d
---

## Solution of a Poisson problem  with discontinuous coefficient in a cube with a sphere inside for the discontinuity

The surface mesh of the cube is built by glueing the 6 meshes of the faces. The mesh of the sphere is built by mapping the mesh of a square

So first the 6 surface meshes for the cube
~~~freefem
verbosity=1;

load "tetgen"
load "medit"
 
// 
meshS ThHex;
real volumetet;  // use in tetg.
{
	//  first  build the 6 faces of the hex.
real x0=-1,x1=1;
real y0=-1.1,y1=1.1;
real z0=-1.2,z1=1.2;

int nx=19,ny=20,nz=21;
//  a  volume  of  on tet. 
volumetet= (x1-x0)*(y1-y0)*(z1-z0)/ (nx*ny*ny) /6.;

mesh Thx = square(ny,nz,[y0+(y1-y0)*x,z0+(z1-z0)*y]);
mesh Thy = square(nx,nz,[x0+(x1-x0)*x,z0+(z1-z0)*y]);
mesh Thz = square(nx,ny,[x0+(x1-x0)*x,y0+(y1-y0)*y]);

int[int] refz=[0,5];  //  down
int[int] refZ=[0,6];   //  up
int[int] refy=[0,3];  //  front
int[int] refY=[0,4];   // back
int[int] refx=[0,1];  // left
int[int] refX=[0,2];   // right


meshS Thx0 = movemesh23(Thx,transfo=[x0,x,y],orientation=-1,region=refx,removeduplicate=false); 
meshS Thx1 = movemesh23(Thx,transfo=[x1,x,y],orientation=1,region=refX,removeduplicate=false);
meshS Thy0 = movemesh23(Thy,transfo=[x,y0,y],orientation=+1,region=refy,removeduplicate=false);
meshS Thy1 = movemesh23(Thy,transfo=[x,y1,y],orientation=-1,region=refY,removeduplicate=false);
meshS Thz0 = movemesh23(Thz,transfo=[x,y,z0],orientation=-1,region=refz,removeduplicate=false);
meshS Thz1 = movemesh23(Thz,transfo=[x,y,z1],orientation=+1,region=refZ,removeduplicate=false);

//medit("  --- ", Thx0,Thx1,Thy0,Thy1,Thz0,Thz1);
 ThHex = Thx0+Thx1+Thy0+Thy1+Thz0+Thz1;
 
}
meshS Thsph;
~~~
Then the surface mesh of the sphere
~~~freefem
{
mesh  Th=square(10,20,[x*pi-pi/2,2*y*pi]);  //  $]\frac{-pi}{2},frac{-pi}{2}[\times]0,2\pi[ $
//  a paratrization of a sphere 
func f1 =cos(x)*cos(y);
func f2 =cos(x)*sin(y);
func f3 = sin(x);
//  de  partiel derivatrive of the parametrization DF
func f1x=sin(x)*cos(y);   
func f1y=-cos(x)*sin(y);
func f2x=-sin(x)*sin(y);
func f2y=cos(x)*cos(y);
func f3x=cos(x);
func f3y=0;
// $  M = DF^t DF $
func m11=f1x^2+f2x^2+f3x^2;
func m21=f1x*f1y+f2x*f2y+f3x*f3y;
func m22=f1y^2+f2y^2+f3y^2;

func perio=[[4,y],[2,y],[1,x],[3,x]];  // to store the periodic condition 

// the intial mesh
//savemesh(Th,"sphere",[f1,f2,f3]);

real R=0.5,hh=0.1/R;// hh  taille du maille sur la shere unite. 
real vv= 1/square(hh);
verbosity=2;
Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,inquire=1,periodic=perio);
plot(Th,wait=1);
Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
plot(Th,wait=1);
Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
plot(Th,wait=1);
Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);

Thsph = movemesh23(Th,transfo=[f1*R,f2*R,f3*R],orientation=-1,removeduplicate=false);
}
~~~
Then define the final volumic mesh and the PDE
~~~freefem
meshS ThS = ThHex+Thsph; // "glueing" all the surface meshes
//medit("Bounday mesh",ThS,wait=1);
// build a mesh of a axis parallel box with TetGen
real[int] domaine = [0,0,0,1,volumetet,0,0,0.7,2,volumetet];
mesh3 Th = tetg(ThS,switch="pqaAAYYQ",nbofregions=2,regionlist=domaine);    
// Tetrahelize the interior of the cube with tetgen
medit("tetg",Th,wait=1);
//savemesh(Th,"Th-hex-sph.mesh");

fespace Ph(Th,P03d);
Ph reg=region;

cout << "  centre = " << reg(0,0,0) << endl;
cout << " exterieur = " << reg(0,0,0.7) << endl;

macro Grad(u) [dx(u),dy(u),dz(u)] // EOM

fespace Vh(Th,P13d);
Vh uh,vh;
real f=1.;
real gn = 1.;
real cf= 1;
problem P(uh,vh)=
   int3d(Th,1)( Grad(uh)'*Grad(vh)*100) 
  +  int3d(Th,2)( Grad(uh)'*Grad(vh)*2) 
  + int3d(Th) (vh*f)
  + on(-1,uh=-1) + on(1,uh=1) 
  + int2d(Th,2,-2)(vh*gn)
  + int2d(Th,3,-3)(cf*vh*uh)
  ; 
  P;
  
// FFCS: with 3D view parameters
real[int] CameraPositionValue = [3.50634,-2.51489,2.60313];
real[int] CameraFocalPointValue = [0.0604689,-0.304636,-0.256484];
real[int] CameraViewUpValue = [0.7198,0.502367,-0.479078];
real[int] CutPlaneOriginValue = [-0.5,-0.55,0.0335184];
real[int] CutPlaneNormalValue = [0,0,1];

plot(uh,wait=1, nbiso=6,
	BorderAsMesh = 1,
	CameraPosition=CameraPositionValue,
	CameraFocalPoint=CameraFocalPointValue,
	CameraViewUp=CameraViewUpValue,
	CutPlaneOrigin=CutPlaneOriginValue,
	CutPlaneNormal = CutPlaneNormalValue);
medit("   uh ",Th, uh,wait=1); 

~~~

| The solution displayed with medit and plot |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/Poisson-cube-ballon/solution.png