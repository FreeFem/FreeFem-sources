load "medit"
int nn=3;
mesh Th2=square(nn,nn);
fespace Vh2(Th2,P2);  Vh2 ux,uz,p2;
int[int] rup=[0,2],  rdown=[0,1], rmid=[1,1,2,1,3,1,4,1];
real zmin=0,zmax=1;
mesh3 Th=buildlayers(Th2,nn,
  zbound=[zmin,zmax],  labelmid=rmid, 
  reffaceup = rup,     reffacelow = rdown);
  
medit("c10x10x10",Th,wait=1);

// FFCS: testing 3d plots
plot(Th);

fespace VVh(Th,[P2,P2,P2,P1]);
fespace Vh(Th,P23d);
macro Grad(u) [dx(u),dy(u),dz(u)]// EOM
macro div(u1,u2,u3) (dx(u1)+dy(u2)+dz(u3)) //EOM

VVh [u1,u2,u3,p];
VVh [v1,v2,v3,q];
func fup = (1-x)*(x)*y*(1-y)*16;
solve vStokes([u1,u2,u3,p],[v1,v2,v3,q]) = 
  int3d(Th,qforder=3)( Grad(u1)'*Grad(v1) +  Grad(u2)'*Grad(v2) +  Grad(u3)'*Grad(v3) //)';
		       - div(u1,u2,u3)*q - div(v1,v2,v3)*p + 1e-10*q*p ) 
  + on(2,u1=fup,u2=0,u3=0) + on(1,u1=0,u2=0,u3=0) ;
plot(p,wait=1, nbiso=5); // a 3d plot of iso  pressure. 
plot([u1,u2,u3] ,wait=1, nbiso=5); // a 3d plot of iso  pressure. 
//  to see the 10 cup plan in 2d 
for(int i=1;i<10;i++)
  {
    real yy=i/10.;
    ux= u1(x,yy,y);
    uz= u3(x,yy,y);
    p2= p(x,yy,y);
    plot([ux,uz],p2,cmm=" cut y = "+yy,wait= 1);
  }

// FFCS: testing 3d plots
plot(u1);
plot(u2);
plot(u3);
plot(p);
