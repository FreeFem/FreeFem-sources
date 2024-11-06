load "medit"// buildlayer

border C(t=0,2*pi) { x = cos(t); y=sin(t); label=1;}
mesh Baseh = buildmesh(C(20));
plot(Baseh,wait=1);

int[int] rup=[0,1],  rdown=[0,2], rmid=[1,3];
func zmin= 1;
func zmax= 10;
int nlayer=100;
mesh3 Th=buildlayers(Baseh,nlayer,
  coef= 1.,
  zbound=[zmin,zmax],
  labelmid=rmid, 
  reffaceup = rup,
  reffacelow = rdown);
medit("Cyl",Th,wait=1);
// FFCS: testing 3d plots
plot(Th,cmm="Cyl");
