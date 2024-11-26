---
name: cone
category: mesh
folder: 3d
---

## Build a volumic mesh of a cone

| Easiest is build a 2D generator curve, here it is a triangle,|
|------------------------|
|![][_Triangle]          |

|and then use a 2D mesh.|
|-----------------------|
|![][_mesh2D]           |

~~~freefem
include "ball-buildlayer.idp"
load "medit"
// cone using buildlayers with a triangle 
real RR=1,HH=1; 
border Taxe(t=0,HH){x=t;y=0;label=0;};
border Hypo(t=1,0){x=HH*t;y=RR*t;label=1;};
border Vert(t=0,RR){x=HH;y=t;label=2;};

int nn=10;
real h= 1./nn;
plot( Taxe(HH*nn)+ Hypo(sqrt(HH*HH+RR*RR)*nn) + Vert(RR*nn), dim=2, wait=1 );
mesh Th2=buildmesh(  Taxe(HH*nn)+ Hypo(sqrt(HH*HH+RR*RR)*nn) + Vert(RR*nn) ) ;
plot(Th2,wait=1);
~~~
 and lift it 3D by  rotation around an axis
~~~freefem
mesh3 Th3T=BuildAxiOx(Th2,h);
//medit("cone",Th3T,wait=1);
plot(Th3T,cmm="cone");
~~~

| The final result: |
|-------------------|
|![][_cone]         |


[_Triangle]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/cone/Triangle.png

[_mesh2D]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/cone/mesh2D.png

[_cone]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/cone/cone.png
