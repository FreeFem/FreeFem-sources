load "medit"
int nn=10;
border cc(t=0,2*pi){x=cos(t);y=sin(t);label=1;}
mesh Th2= buildmesh(cc(100));
fespace Vh2(Th2,P2);
Vh2 ux,uz,p2;
int[int] rup=[0,2],  rdown=[0,1], rmid=[1,1,2,1,3,1,4,1];
func zmin= 2-sqrt(4-(x*x+y*y));
func zmax= 2-sqrt(3.);

mesh3 Th=buildlayers(Th2,nn,
  coef=  (zmax-zmin)/zmax ,
  zbound=[zmin,zmax],
  labelmid=rmid, 
  reffaceup = rup,
  reffacelow = rdown);
medit("Lac",Th,wait=1);
// FFCS: testing 3d plots
plot(Th,cmm="Lac");
