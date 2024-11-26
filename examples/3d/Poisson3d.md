// build de mesh of a Sphere
// -------------------------- 
load "tetgen"
load "medit"

mesh Th=square(10,20,[x*pi-pi/2,2*y*pi]);  //  $]\frac{-pi}{2},frac{-pi}{2}[\times]0,2\pi[ $
//  a parametrization of a sphere 
func f1 =cos(x)*cos(y);
func f2 =cos(x)*sin(y);
func f3 = sin(x);
//  de  partiel derivative of the parametrization DF
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

func perio=[[4,y],[2,y],[1,x],[3,x]];  
real hh=0.1;
real vv= 1/square(hh);
verbosity=2;
Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
plot(Th,wait=1);

verbosity=2;
real[int] domaine =[0.,0.,0.,1,0.01];
mesh3 Th3=tetgtransfo(Th,transfo=[f1,f2,f3],nbofregions=1,regionlist=domaine);
//savemesh(Th3,"sphere.meshb");
medit("sphere",Th3);

// FFCS - check 3D plots
plot(Th3,cmm="sphere");

fespace Vh(Th3,P23d);
func ue =   2*x*x + 3*y*y + 4*z*z+ 5*x*y+6*x*z+1;
func f= -18. ;
Vh uhe = ue; // bug ..
cout << " uhe min:  " << uhe[].min << " max:" << uhe[].max << endl;
cout << uhe(0.,0.,0.) << endl;


//savesol("f3.sol",Th3,ue,ue,[f,ue,f],order=1);
//int bb=meditmeshsol("sol",Th3,solution=1,scalar=uhe);


border cc(t=0,2*pi){x=cos(t);y=sin(t);label=1;}
mesh Th2=buildmesh(cc(50));
fespace Vh2(Th2,P2);

Vh u,v;

macro Grad3(u) [dx(u),dy(u),dz(u)]  // EOM

problem Lap3d(u,v,solver=CG)=int3d(Th3)(Grad3(v)' *Grad3(u)) - int3d(Th3)(f*v) + on(0,1,u=ue);
Lap3d;
cout << " u min::   " << u[]. min << "  max: " << u[].max << endl;
real err= int3d(Th3)( square(u-ue) );
cout << int3d(Th3)(1.) << " = " << Th3.measure << endl;
Vh d= ue-u;
cout <<  " err = " << err <<  " diff l^\intfy = " << d[].linfty << endl;
Vh2 u2=u,u2e=ue;
plot(u2,wait=1);
plot(u2,u2e,wait=1);

// FFCS - check 3D plots
plot(u);

assert(err < 1e-9);
