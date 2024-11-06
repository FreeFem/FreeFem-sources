//  example to build a mesh a cone 
load "medit"
// cone using buildlayers with a triangle 
real LX=1,LY=1,LXY=sqrt(LX*LX+LY*LY),HH=1; 
border Hypo(t=1,0){x=LX*t;y=LY*(1-t);label=1;};
border Vert(t=LY,0){x=0;y=t;label=0;};
border Hori(t=0,LX){x=t;y=0;label=0;};

int nn=10;
real h= 1./nn;
plot(Vert(LY*nn)+ Hypo(LXY*nn) + Hori(LX*nn),wait=1);
mesh Th2=buildmesh( Vert(LY*nn)+ Hypo(LXY*nn) + Hori(LX*nn) ) ;
 Th2 = Th2 + movemesh(Th2,[x,-y])+ movemesh(Th2,[-x,-y])+  movemesh(Th2,[-x,y]);
plot(Th2,wait=1);
func fpyramide= (1-abs(x)/LX-abs(y)/LY)*HH;
fespace Vh2(Th2,P1);
Vh2 fp2=fpyramide;
plot(fp2,wait=1,dim=3);


int[int] r1T=[0,0], r2up=[0,1],r2down=[0,1];
int[int] r4T=[0,2]; 
mesh3 Th3=buildlayers(Th2,coef= max(fpyramide/HH,0.01), nn,zbound=[0,fpyramide],
 region=r1T, labelup=r2up, labeldown=r2down);

medit("Pyramide",Th3,wait=1);
// FFCS: testing 3d plots
plot(Th3);
