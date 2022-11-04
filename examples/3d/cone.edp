//  example to build a mesh a cone 
include "ball-buildlayer.idp"
load "medit"
// cone using buildlayers with a triangle 
real RR=1,HH=1; 
border Taxe(t=0,HH){x=t;y=0;label=0;};
border Hypo(t=1,0){x=HH*t;y=RR*t;label=1;};
border Vert(t=0,RR){x=HH;y=t;label=2;};

int nn=10;
real h= 1./nn;
mesh Th2=buildmesh(  Taxe(HH*nn)+ Hypo(sqrt(HH*HH+RR*RR)*nn) + Vert(RR*nn) ) ;

mesh3 Th3T=BuildAxiOx(Th2,h);
medit("cone",Th3T,wait=1);
// FFCS: testing 3d plots
plot(Th3T,cmm="cone");
