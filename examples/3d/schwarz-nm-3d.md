---
name: schwarz-nm-3d
category: math
layout: 3d
---

## Use metis to split a triangulation then solve a Laplacian by Schwarz DDM

The mesh of a cube is split using metis.
~~~freefem
bool withmetis=1;
bool RAS=0;
int sizeoverlaps=2; // size off overlap 
int nnx=2,nny=2,nnz=2;

func bool AddLayers(mesh3 & Th,real[int] &ssd,int n,real[int] &unssd)
{
  //  build a continuous function  uussd (P1) :
  //  ssd in the caracteristics function on the input sub domain.
  //  such that : 
  //   unssd = 1 when   ssd =1;
  //   add n layer of element (size of the overlap)
  //   and unssd = 0 ouside of this layer ...
  // ---------------------------------
  fespace Vh(Th,P1);
  fespace Ph(Th,P0);
  Ph s;
  assert(ssd.n==Ph.ndof);
  assert(unssd.n==Vh.ndof);
  unssd=0;
  s[]= ssd;
  //  plot(s,wait=1,fill=1);
  Vh u;
  varf vM(u,v)=int3d(Th,qforder=1)(u*v/volume);
  matrix M=vM(Ph,Vh);
  
  for(int i=0;i<n;++i)
    {
      u[]= M*s[];
      // plot(u,wait=1);
      u = u>.1; 
      // plot(u,wait=1);
      unssd+= u[];
      s[]= M'*u[];//'
      s = s >0.1;
    }
  unssd /= (n);
  u[]=unssd;
  ssd=s[];      
  return true;
}

int withplot=3;
include "cube.idp" 
 int[int]  NN=[25,25,25]; //  the number of step in each  direction                                                                                                                    
 real [int,int]  BB=[[0,1],[0,1],[0,1]]; // bounding box                                                                                                             
 int [int,int]  L=[[1,1],[1,1],[1,1]]; // the label of the 6 face left,right,
//  front, back, down, right
mesh3 Th=Cube(NN,BB,L);
int npart= nnx*nny*nnz;
fespace Ph(Th,P0);
fespace Vh(Th,P1);

Ph  part;
Vh  sun=0,unssd=0;
Ph xx=x,yy=y,zz=z,nupp;
//part = int(xx*nnx)*nny + int(yy*nny);
part = int(xx*nnx)*nny*nnz + int(yy*nny)*nnz + int(zz*nnz);
//plot(part,wait=1);
if(withmetis)
  {
    load "metis";
    int[int] nupart(Th.nt);
    metisdual(nupart,Th,npart); 
    for(int i=0;i<nupart.n;++i)
      part[][i]=nupart[i];
  }
if(withplot>1)
plot(part,fill=1,cmm="dual",wait=1);
~~~
Then the PDEs are prepared
~~~freefem
mesh3[int] aTh(npart);
mesh3 Thi=Th;
fespace Vhi(Thi,P1);
Vhi[int] au(npart),pun(npart);
matrix[int] Rih(npart);
matrix[int] Dih(npart);
matrix[int] aA(npart);
Vhi[int] auntgv(npart);
Vhi[int] rhsi(npart);

for(int i=0;i<npart;++i)
  {
    Ph suppi= abs(part-i)<0.1;
    AddLayers(Th,suppi[],sizeoverlaps,unssd[]);
    Thi=aTh[i]=trunc(Th,suppi>0,label=10,split=1);
    Rih[i]=interpolate(Vhi,Vh,inside=1); //  Vh -> Vhi
    if(RAS)
      {
        suppi= abs(part-i)<0.1;
        varf vSuppi(u,v)=int3d(Th,qforder=1)(suppi*v/volume);
        unssd[]= vSuppi(0,Vh);
        unssd = unssd>0.;
        if(withplot>19)
          plot(unssd,wait=1);
      }
    pun[i][]=Rih[i]*unssd[];
    sun[] += Rih[i]'*pun[i][];//';
    if(withplot>9)
      plot(part,aTh[i],fill=1,wait=1);
  }
real[int] viso=[0,0.1,0.2,0.3];  
plot(sun,wait=1,dim=3,fill=1,viso=viso);
for(int i=0;i<npart;++i)
  {
    Thi=aTh[i];
    pun[i]= pun[i]/sun;
    if(withplot>8)
      plot(pun[i],wait=1);    
  }

//  verif partition of unite 

macro Grad(u) [dx(u),dy(u),dz(u)]//EOM 
  sun=0;

for(int i=0;i<npart;++i)
  {
    cout << " build part :" << i << "/" << npart << endl;
    Thi=aTh[i];
    varf va(u,v) = 
      int3d(Thi)(Grad(u)'*Grad(v))//')
      +on(1,u=1) + int3d(Thi)(v)
      +on(10,u=0) ; 
    
    cout << i << " -----------Vhi.ndof " << Vhi.ndof << endl;
    aA[i]=va(Vhi,Vhi);
      cout << i << " -----------Vhi.ndof " << Vhi.ndof << endl;
  
    set(aA[i],solver="SPARSESOLVER");
    rhsi[i][]= va(0,Vhi);
    Dih[i]=pun[i][];
    real[int]  un(Vhi.ndof);
    un=1.;
    real[int] ui=Dih[i]*un; 
    sun[] += Rih[i]'*ui;;//';
    varf vaun(u,v) = on(10,u=1);
    auntgv[i][]=vaun(0,Vhi); // store array of tgv on Gamma intern.
  }
if(withplot>5)
  plot(sun,fill=1,wait=1);
cout << sun[].max << " " << sun[].min<< endl;
// verification of the partition of the unite.
assert( 1.-1e-9 <= sun[].min  && 1.+1e-9 >= sun[].max);  
~~~
The Schwarz iterations are implemented
~~~freefem
int nitermax=1000;
{
  Vh un=0;
  for(int iter=0;iter<nitermax;++iter)
    {
      real err=0;
      Vh un1=0;
      for(int i=0;i<npart;++i)
        {
          Thi=aTh[i];
          real[int] ui=Rih[i]*un[];//';
          //{   Vhi uuu; uuu[]=ui;      plot(uuu,wait=1);}
          real[int] bi = ui .* auntgv[i][];
          bi = auntgv[i][] ? bi :  rhsi[i][];  
          ui=au[i][];
          ui= aA[i] ^-1 * bi;
          //{   Vhi uuu; uuu[]=ui;      plot(uuu,wait=1);}
          bi = ui-au[i][];
          err += bi'*bi;//';
          au[i][]= ui;
          bi = Dih[i]*ui;
          un1[] += Rih[i]'*bi;//';
        }
      err= sqrt(err);
      cout << iter << " Err = " << err << endl;
      if(err<1e-3) break;
      //    plot(un1,wait=1);
      un[]=un1[];
      if(withplot>2)
        plot(au,dim=3,wait=0,cmm=" iter  "+iter,fill=1 );
    }
  plot(un,wait=1,dim=3,fill=1);
}
~~~

|The solution            |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/schwarz-nm-3d/solution.png