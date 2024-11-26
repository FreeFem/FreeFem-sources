// Sphere from a  partir d'un icosahedron
// Guillaume Vergez. 
load "medit" 
load "tetgen"
load "mmg"


include "MeshSurface.idp"
real R= 2;
real hsize= 0.1; 
real areaT = hsize*hsize*sqrt(3)/2.;
real areaS= 4*pi*R*R;
real nn = sqrt(areaS/areaT/20.); 
cout << nn << endl; 
meshS Thb=Sphere20(2.,nn,1,0);// sphere of  Raduis 2 ..
real[int] bb(6);
boundingbox(Thb,bb);
cout << " BB = " << bb << endl;
plot(Thb,wait=1);//Thb=mmgs(Thb,hmin=hsize,hmax=hsize,hgrad=1.1);

real[int] domain = [0.,0.,0.,1,hsize^3/6];
mesh3 Th3sph=tetg(Thb,switch="paAAQYY",nbofregions=1,regionlist=domain);


 plot(cmm="sphere",Th3sph,wait=1);
