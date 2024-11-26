load "medit"
int nn=10;
// a pententiel flow on a lac .. 
//  a first freefem++ 3d example 
// ------  not to bad ......
verbosity=3;
mesh Th2("lac-leman-v4.msh");
fespace Vh2(Th2,P1);
Vh2 deep;
{  Vh2 v; 
	macro Grad(u) [dx(u),dy(u)] //
	solve P(deep,v)= int2d(Th2)(Grad(deep)'*Grad(v))+int2d(Th2)(v)
	+on(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,deep=-1);
	deep = deep*5/abs(deep[].min);
	plot(deep,wait=1);
}
Vh2 ux,uz,p2;
int[int] rup=[0,200],  rdown=[0,100],rmid(17*2);
for(int i=0;i<rmid.n;++i)
  rmid[i]=1+i/2;
cout << rmid << endl;
real maxdeep = deep[].min;
mesh3 Th=buildlayers(Th2,nn,
  coef= deep/maxdeep,
  zbound=[deep,0],
  labelmid=rmid, 
  reffaceup = rup,
  reffacelow = rdown);
if(!NoUseOfWait)medit("Leman",Th);
fespace Vh(Th,P13d);
Vh p,q;
//  (-deep[].min)*c = 0.5 Km
//  c = 
real cc=(0.5/-deep[].min);
cout << cc << " cc = " << endl;
cc=1; //  otherwise bug in boundaing condition ... 
macro Grad(u) [dx(u),dy(u),cc*dz(u)] //

real ain=int2d(Th,1)(1.);
real aout=int2d(Th,2)(1.);
cout << " area " << ain << " " << aout << endl;
real din=1./ain;
real dout=-1./aout;
solve P(p,q)= int3d(Th)(Grad(p)'*Grad(q)+1e-5*p*q)-int2d(Th,1)(q*din)+int2d(Th,2)(q*dout);

plot(p,wait=1,nbiso=30,value=1);
if(!NoUseOfWait) medit("potentiel",Th,p,wait=1);

real vregtest = p[].max ;
cout << " p max " << vregtest << endl; 
