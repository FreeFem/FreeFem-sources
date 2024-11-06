/*
   $ - \Delta p = f $   on $\Omega$, 
   $ dp / dn = (g1d,g2d). n  $ on $\Gamma_{1}$ 
   $ p = gd  $ on $\Gamma_{2}$    
   with de Mixte finite element formulation 

   Find $p\in L^2(\Omega) $  and $u\in H(div) $ such than 
   $$  u - Grad p = 0    $$
   $$ - div u =  f $$
   $$  u. n = (g1d,g2d). n   \mbox{ on } \Gamma_{2}$$
   $$ p = gd  \mbox{ on }\Gamma_{1}$$
   the variationnel form is: 
                                                                                                                   
 $\forall v\in H(div)$;  $v.n = 0$ on $\Gamma_{2} $:  

  $ \int_\Omega  u v + p div v -\int_{\Gamma_{1}} gd* v.n  = 0 $ 
 $\forall q\in L^2$:   $  +\int_\Omega q div u = -\int_Omega f q  $
and $ u.n = (g1n,g2n).n$ on $\Gamma_2$

*/
include "cube.idp"
    int[int]  Nxyz=[10,10,10];
    real [int,int]  Bxyz=[[0,1],[0,1],[0,1]];
    int [int,int]  Lxyz=[[1,1],[1,1],[2,1]];
mesh3 Th=Cube(Nxyz,Bxyz,Lxyz);
fespace Vh(Th,P1);
fespace Rh(Th,RT03d);
fespace Nh(Th,Edge03d);//  Nedelec Finite element. 
fespace Ph(Th,P0);

func gd = 1.;

func g1n = 2.;
func g2n = 3.; 
func g3n = 4.; 

func f = 1.;

Rh [u1,u2,u3],[v1,v2,v3];
Nh [e1,e2,e3];
[u1,u2,u3]=[1+100*x,2+100*y,3+100*z];

// a + b ^ x = 
/*
  b1    x     a1 + b2*z - b3*y 
  b2 ^  y  =  a2 - b1*z + b3*x
  b3    z     a3 + b1*y - b2*x
*/
real b1=30,b2=10,b3=20;
func ex1=100+b2*z-b3*y;

func ex1x=0.;
func ex1y=-b3+0;
func ex1z=b2+0;

func ex2=200.- b1*z + b3*x ;
func ex2x= b3 +0;
func ex2y= 0. ;
func ex2z= -b1 +0;
func ex3=300.+b1*y - b2*x ;
func ex3x= -b2 +0;
func ex3y= b1 +0;
func ex3z= 0. ;
[e1,e2,e3]=[ex1,ex2,ex3]; 

int k=Th(.1,.2,.3).nuTriangle ;
cout << " u = " << u1(.1,.2,.3)  << " " << u2(.1,.2,.3) << " " << u3(.1,.2,.3) << endl;
cout << " dx u = " << dx(u1)(.1,.2,.3)  << " " << dy(u2)(.1,.2,.3) << " " << dz(u3)(.1,.2,.3) << endl;

cout << " e  = " << e1(.1,.2,.3)  << " " << e2(.1,.2,.3) << " " << e3(.1,.2,.3) << endl;
cout << " ex = " << ex1(.1,.2,.3)  << " " << ex2(.1,.2,.3) << " " << ex3(.1,.2,.3) << endl;


cout << " dx,dy,dz   e1x= " << ex1x(.1,.2,.3)  << " " << ex1y(.1,.2,.3) << " " << ex1z(.1,.2,.3) << endl;
cout << " dx,dy,dz   e2x= " << ex2x(.1,.2,.3)  << " " << ex2y(.1,.2,.3) << " " << ex2z(.1,.2,.3) << endl;
cout << " dx,dy,dz   e3x= " << ex3x(.1,.2,.3)  << " " << ex3y(.1,.2,.3) << " " << ex3z(.1,.2,.3) << endl;

cout << " dx,dy,dz   e1 = " << dx(e1)(.1,.2,.3)  << " " << dy(e1)(.1,.2,.3) << " " << dz(e1)(.1,.2,.3) << endl;
cout << " dx,dy,dz   e2 = " << dx(e2)(.1,.2,.3)  << " " << dy(e2)(.1,.2,.3) << " " << dz(e2)(.1,.2,.3) << endl;
cout << " dx,dy,dz   e3 = " << dx(e3)(.1,.2,.3)  << " " << dy(e3)(.1,.2,.3) << " " << dz(e3)(.1,.2,.3) << endl;


cout << " k = " << k << endl;
cout << Rh(k,0) << " " <<Rh(k,1) << " " <<Rh(k,2) << " " <<Rh(k,3) << endl;
cout << " df = " << u1[][Rh(k,0)] <<  " " << u1[][Rh(k,1)]  <<" " << u1[][Rh(k,2)]  << " " << u1[][Rh(k,2)] << endl;
// cout << u1[] << endl;

Vh P,Q;
Ph p,q; 
macro div(u1,u2,u3) (dx(u1)+dy(u2)+dz(u3)) //
macro Grad(u) [dx(u),dy(u),dz(u)]  //
  problem laplace(P,Q,solver=CG) = 
  int3d(Th) ( Grad(P)'*Grad(Q)) //') for emacs
  - int3d(Th)(f*Q) 
  + on(1,P=gd) 
  - int2d(Th,2) ( (g1n*N.x+g2n*N.y+g3n*N.z)*Q);

fespace RPh(Th,[RT03d,P0]);
varf von1([u1,u2,u3,p],[v1,v2,v3,q])  = 
   int3d(Th)( p*q*1e-15+ u1*v1 + u2*v2 + u3*v3 + p*div(v1,v2,v3) + div(u1,u2,u3)*q )
 - int3d(Th) ( f*q)
 + int2d(Th,1)( gd*(v1*N.x +v2*N.y + v3*N.z) )  //  int on gamma 
 + on(2,u1=g1n,u2=g2n,u3=g3n);

RPh [vv1,vv2,vv3,qq];
// some verification Boundary Condition
// and interpolation ...
real[int]  ron=von1(0,RPh);
vv3[]=von1(0,RPh);
cout << " vv: = " << vv1(.1,.2,.001)  << " " << vv2(.1,.2,.001) << " " << vv3(.1,.2,.001) << endl;
[vv1,vv2,vv3,qq]=[g1n,g2n,g3n,100];
[v1,v2,v3]=[g1n,g2n,g3n];

cout << " vv: = " << vv1(.1,.2,.001)  << " " << vv2(.1,.2,.001) << " " << vv3(.1,.2,.001) << " " << qq(.1,.2,.001) << endl;
cout << " v : = " << v1(.1,.2,.001)  << " " << v2(.1,.2,.001) << " " << v3(.1,.2,.001)  << endl;

// end of verification of Boundary Condition ... 

problem laplaceMixte([u1,u2,u3,p],[v1,v2,v3,q],eps=1.0e-10,tgv=1e30,dimKrylov=1000) =
   int3d(Th)( p*q*1e-15+ u1*v1 + u2*v2 + u3*v3 + p*div(v1,v2,v3) + div(u1,u2,u3)*q )
 + int3d(Th) ( f*q)
 - int2d(Th,1)( gd*(v1*N.x +v2*N.y + v3*N.z) )  //  int on gamma 
 + on(2,u1=g1n,u2=g2n,u3=g3n);

laplace;

// FFCS: add 3D view parameters
real[int] CameraPositionValue = [0.0165449,3.23891,-0.991528];
real[int] CameraFocalPointValue = [0.5,0.5,0.5];
real[int] CameraViewUpValue = [0.671735,0.442219,0.594318];
real[int] CutPlaneOriginValue = [0.5,0.5,1.01];
real[int] CutPlaneNormalValue = [0.689523,0.722423,0.0516115];
plot(P,fill=0,boundary=0,ShowMeshes=1,CutPlane=1,
	CameraPosition=CameraPositionValue,
	CameraFocalPoint=CameraFocalPointValue,
	CameraViewUp=CameraViewUpValue,
	CutPlaneOrigin=CutPlaneOriginValue,
	CutPlaneNormal = CutPlaneNormalValue);

laplaceMixte;

real errL2=sqrt(int3d(Th)(square(P-p))) ;
cout << " int 2 x,yz "<<int2d(Th,2)(x) << " " << int2d(Th,2)(y) << " " << int2d(Th,2)(z) << endl;
cout << " int 2 gn "<<int2d(Th,2)(g1n) << " " << int2d(Th,2)(g2n) << " " << int2d(Th,2)(g3n) << endl;
cout << " int 2 U  "<<int2d(Th,2)(u1) << " " << int2d(Th,2)(u2) << " " << int2d(Th,2)(u3) << endl;
cout << " int 2 V  "<<int2d(Th,2)(vv1) << " " << int2d(Th,2)(vv2) << " " << int2d(Th,2)(vv3) << endl;
cout << " int 2 DP "<<int2d(Th,2)(dx(P)) << " " << int2d(Th,2)(dy(P)) << " " << int2d(Th,2)(dz(P)) << endl;
  
cout << "  diff: u Gamma_2 " <<    sqrt(int2d(Th,2) ( square((g1n*N.x+g2n*N.y+g3n*N.z) - (u1*N.x +u2*N.y + u3*N.z) ) ) ) <<endl;
cout << "  diff: P Gamma_2 " <<    sqrt(int2d(Th,2) ( square((g1n*N.x+g2n*N.y+g3n*N.z) - (dx(P)*N.x +dy(P)*N.y + dz(P)*N.z) ) ) ) <<endl;
cout << " diff err L2 :" << errL2 << endl;
cout << "    P     L2 :" <<sqrt(int3d(Th)(square(P))) << endl;
cout << "    p     L2 :" <<sqrt(int3d(Th)(square(p))) << endl;
assert(errL2<0.05);
