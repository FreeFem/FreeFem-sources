include "ball-buildlayer.idp"

load "medit"
{
real RR=4.8828; 
real h= 0.2; // mesh size 
mesh3 Th=BuildBall(RR,h,2);
cout << " region = "<< regions(Th)<< endl;
cout << " label = "<< labels(Th)<< endl;
medit("Th",Th);
plot(Th,wait=1);
}
{
real RR=4.8828; 
real h= 0.2; // mesh size 

border Taxe(t=-RR,RR){x=t;y=0;label=0;}
border CC(t=0,pi){x=RR*cos(t);y=RR*sin(t);label=2+(x>0);}
mesh Th2=buildmesh(  Taxe(2*RR/h)+ CC(pi*RR/h) ) ;
mesh3 Th=BuildAxiOx(Th2,h);
cout << " region = "<< regions(Th)<< endl;
cout << " label = "<< labels(Th)<< endl;
medit("Th",Th);
plot(Th,wait=1);
}
