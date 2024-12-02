---
name: bottle
category: Solid Mechanics
folder: 3d
author: F. Hecht
---
## Build the mesh of a bottle
The bottle is axisymmetric and has a neck similar to those of wine bottles. The program constructs a closed (the bottle has a cork)  surface mesh Ths. First one construct the generatrix using spline interpolations of polygonal curves from the library GSL (Gnu Scientific Library) and then rotate it around the z-axis.

| The end product: a bottle |
|---------------------------|
|![][_bottle]               |

~~~freefem
load "medit"
load "gsl"
load "tetgen"
real hh=.15;// mesh size

real Htube=8; 
real Hneck=2;
real Hbot=Htube+Hneck;

real Rext = 1.;// External radius
real Rneck= 0.2; // Radius of the neck of the bootle
int labtop =3, labbottom =2, labcyl = 1;

meshS Ths;
{
meshS Th3c,Th3bottom, Th3top;
{
// form of the neck ... 	
real[int,int] srneck= 
[
 [Htube-0.001,Htube,Htube+Hneck*0.1, Htube+Hneck*0.3, Htube+Hneck*0.7, Htube+Hneck*0.9, Htube+Hneck+0.1], [Rext,Rext ,Rext, Rext*.7+Rneck*.3 , Rext*.1 + Rneck*.9, Rneck, Rneck  ]
 ];
gslspline rneck(gslinterpcspline,srneck);
~~~
 rneck is the curve of the bollte form z = 0 to Hbot neck of the bottle
  One needs to adapt the 2D mesh for a better 3D mesh. Th2c  will be the surface minus top and bottom.
~~~freefem
mesh Th2c= square(Hbot/hh,2*pi*Rext/hh,[x*Hbot,y*2*pi]);
fespace V2x(Th2c,P1);
macro Dx(u) [dx(u#1),dx(u#2),dx(u#3)] // u#1 is the concatenated value of u with 1
macro Dy(u) [dy(u#1),dy(u#2),dy(u#3)] //
// The transformation 
func E1 = rneck(x)*cos(y);
func E2 = rneck(x)*sin(y);
func E3 = x;
V2x ex1,ex2,ex3;
// the metric for adapting the mesh
real hh2 = hh*hh;
func  em11 = Dx(ex)'*Dx(ex)/hh2;
func  em21 = Dx(ex)'*Dy(ex)/hh2;
func  em22 = Dy(ex)'*Dy(ex)/hh2; //'
func perio=[[1,x],[3,x]];
for(int i=0;i<4;++i)
{
 ex1=E1;ex2=E2;ex3=E3; 
 Th2c=adaptmesh(Th2c,em11,em21,em22,IsMetric=1,periodic=perio,nbvx=100000);
}
Th2c=change(Th2c,fregion=labcyl); 
plot(Th2c,wait=1);
~~~

| The mesh Th2c shown below is used to build an axisymmetric body |
|----------------------|
|![][_mesh2d]          |

~~~freefem
Th3c = movemesh23(Th2c,transfo=[E1 , E2 , E3]);// maillage exterieur 
plot(Th3c,wait=1);
//medit("Th3c",Th3c);
}
~~~

| The bottle is open on both ends |
|---------------------------------|
|![][_bottlezero]                 |


The bottom and the top are constructed by gluing 2 meshes of disks but the points on the border must be identical to the borders of the lateral surface
~~~freefem
// extraction of the border of the bottle
int[int] databoder(1);
int ncb= getborder(Th3c,databoder); 
int ktop= Th3c(databoder[databoder[1]]).z < Th3c(databoder[databoder[0]]).z; 
int  kbot=1-ktop;
real Zbot = Th3c(databoder[databoder[kbot]]).z; 
real Ztop = Th3c(databoder[databoder[ktop]]).z; 
cout << " Z  top  " << Zbot  << endl; 
cout << " Z  bot  " << Ztop << endl; 
assert( ncb ==2);
macro DefBorder(bname,kk,Th3,bb,ll)
int n#bname= bb[kk+1]-bb[kk];
border bname(t=bb[kk], bb[kk+1])
{
	real iv = int(t);
	if( iv == bb[kk+1]) iv = bb[kk];

    cout << t << " iv = " << iv << endl;
		iv = bb[iv];
	x= Th3(iv).x ;
	y= Th3(iv).y ;
	label = ll;	
}//
DefBorder(btop,ktop,Th3c,databoder,1)
DefBorder(bbot,kbot,Th3c,databoder,1)
cout << " btop " << nbtop << " bot: " << nbbot << endl;

plot(btop(nbtop)+bbot(-nbbot),wait=1);
Th3bottom=movemesh23(change(buildmesh(bbot(nbbot),fixeborder=1),fregion=labbottom),transfo=[x,y,Zbot],orientation=-1);
Th3top=movemesh23(change(buildmesh(btop(-nbtop),fixeborder=1),fregion=labtop),transfo=[x,y,Ztop],orientation=1);
Ths = Th3c + Th3bottom + Th3top; 
}

real[int] domaine = [0,0,Htube,1,hh^3/6.];
mesh3 Th=tetg(Ths,switch="pqaAYY",regionlist=domaine);
medit("Th",Th);
~~~
Medit can rotate the object. Press "esc" to end medit.


[_bottle]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/bottle/bottle.png


[_mesh2d]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/bottle/mesh2d.png

[_bottlezero]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/bottle/bottlezero.png
