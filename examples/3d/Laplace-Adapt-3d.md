---
name: Laplace-Adapt-3d
category: Applied Math
folder: 3d
---

## Solve the Laplace Equations in a Cube with Discontinous Galerkin Method of degree 1 and mesh refinement

This is a continuation of LapDG3d1.md. The problem is the same
$$
-\Delta u = f,\texttt{ in } \Omega\quad u|_{\partial\Omega}=g
$$
where $\Omega$ is the unit cube minus a half unit cube, $f=1$, $g=1$.
The finite element space chosen if the discontinuous $P^1$.

The unit cube is constructed in a different way with $\texttt{buildlayers}$.

~~~freefem
load "tetgen"
load "mshmet"
load "medit"
//build initial mesh
int nn  = 6;
int[int] l1111=[1,1,1,1],l01=[0,1],l11=[1,1];//   label numbering 
mesh3 Th3=buildlayers(square(nn,nn,region=0,label=l1111),
      nn,  zbound=[0,1],  labelmid=l11,   labelup = l01,  labeldown = l01);
Th3 = trunc(Th3,(x<0.5) | (y < 0.5) | (z < 0.5) ,label=1);// remove the $]0.5,1[^3 cube$
//end of build initial mesh
fespace Vh(Th3,P1);
Vh u,v,usol;

macro Grad(u) [dx(u),dy(u),dz(u)] // EOM

problem Poisson(u,v,solver=CG) = int3d(Th3)( Grad(u)'*Grad(v) )  // ') for emacs 
  -int3d(Th3)( 1*v ) + on(1,u=0);

real errm=1e-2; // level of error

~~~

The  function $\texttt{mshmet}$ refine a mesh from the given Th3 mesh using a metric constructed from the Hessian of $u$, to be used in conjunction with $\texttt{tetgreconstruction}$.
~~~freefem
for(int ii=0; ii<5; ii++)
{
  Poisson;
  cout <<" u min, max = " <<  u[].min << " "<< u[].max << endl;
  Vh h ;
  h[]=mshmet(Th3,u,normalization=1,aniso=0,nbregul=1,hmin=1e-3,hmax=0.3,err=errm);//loptions=MSHloptions,doptions=MSHdoptions);
  cout <<" h min, max = " <<  h[].min << " "<< h[].max << " " << h[].n << " " << Th3.nv << endl;
  // FFCS: add 3D view parameters
  plot(u,wait=1,fill=0,boundary=0,CutPlane=0,ShowMeshes=1,LabelColors=0);
  errm*= 0.8;// change the level of error
  cout << " Th3" << Th3.nv < " " << Th3.nt << endl;
  Th3=tetgreconstruction(Th3,switch="raAQ",sizeofvolume=h*h*h/6.);
  medit("U-adap-iso-"+ii,Th3,u,wait=1);
}
~~~

| The initial mesh       |
|------------------------|
|![][_solution1]         |

|The mesh after 3 refinement |
|----------------------------|
|![][_solution2]             |

| The mesh after  5 refinement |
|------------------------------|
|![][_solution4]               |

| The solution is recomputed after 5 refinement |
|-------------------------|
|![][_solution5]          |

[_solution1]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/Laplace-Adapt-3d/solution1.png

[_solution2]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/Laplace-Adapt-3d/solution2.png

[_solution4]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/Laplace-Adapt-3d/solution4.png

[_solution5]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/Laplace-Adapt-3d/solution5.png
