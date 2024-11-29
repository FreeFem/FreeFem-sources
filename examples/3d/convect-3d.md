---
name: convect-3d
category: fluid
layout: 3d
---

##  Convection of a Tracer
The concentration $v$ of a tracer in a fluid $\Omega$ moving at velocity ${\bf u}=[u_1,u_2,u_3]^T$ is modeled by the convection equation
$$
\partial_t v + {\bf u}\cdot\nabla v=0 \texttt{ in } \Omega\times(0,T)
$$
Here $\Omega$ is a cube, built by extruding a square vertically using $\texttt{buildlayers}$
$v|_{t=0}=v^0$ must be given as well as $v|_\Sigma= v^\Sigma$ where
$$
\Sigma=\{x\in\partial\Omega, t\in(0,T) ~:~{\bf n}_\Gamma\cdot{\bf u} <0\}.
$$

```freefem
int nn=8;

mesh Th2=square(nn,nn,[x*2.-1.,y*2.-1.]);
fespace Vh2(Th2,P1);
int[int] rup=[0,2],  rdown=[0,1], rmid=[1,1,2,1,3,1,4,1];
```
These arrays are used by buildlayers to assign labels to faces. Hence because of rup[0,2] the region 0 in the triangulation of the square becomes the top face of the cube with label = 2.  Similarly the bottom face has label 1 because of rdown[] and finally the sides  which had label 1,2,3,4 (because of the boundary labels of $\texttt{square}$ become faces with label=1.
```freefem
real zmin=-1,zmax=1.;

mesh3 Th=buildlayers(Th2,nn,
  zbound=[zmin,zmax],
  // region=r1, 
  labelmid=rmid, 
  reffaceup = rup,
  reffacelow = rdown);
```
The concentration of the tracer is shaped like a hill around the $[x_0,y_0,z_0]^T$ at time $t=0$.  The velocity is $[1,2,3]^T$.
```freefem
func  real hill(real r2){return exp(-10.*(r2));};

fespace Vh(Th,P13d);

macro Grad(u) [dx(u),dy(u),dz(u)]// EOM
macro div(u1,u2,u3) (dx(u1)+dy(u2)+dz(u3)) //EOM

Vh v,vo;
Vh2 v2;
real x0=0.3,y0=0.3,z0=0;
vo=hill(square(x-x0)+square(y-y0)+square(z-z0));

real t=0;
v2=vo(x,y,0);
plot(v2,cmm=" cut y = 0.5, time ="+t,wait=1);
real dt=0.1;
func u1=1.;
func u2=2.;
func u3=3.;
```
The PDE is discretized by the method of characteristics
$$
v_h(x,t+\delta t)= \Pi v_h(x-\delta t{\bf u}(x,t),t).
$$
where $\Pi$ denotes the $P^1$ interpolation operator.
The method is simple but diffusive. A more precise, elbeit more expensive, discretization is
$$
\int_\Omega v_h(x,t+\delta t) \hat v dx = \int_\Omega v(x-\delta t{\bf u}(x,t),t) \hat v dx ~~\text{ for all hat function of }V_h.
$$
```freefem
verbosity = 1;
v=convect([u1,u2,u3],-dt,vo);
verbosity = 1;
v2=v(x,y,0);
t += dt;
plot(v2,cmm=" cut y = 0.5, time ="+t,wait=1);
// verification ...
int err=0; 
macro Verif(w,val)
{
   real so= int3d(Th)(vo);
   real soi= int3d(Th)(vo*w);
   real sv= int3d(Th)(v);
   real svi= int3d(Th)(v*w);

   cout  << Stringification(w) << "  old " <<  soi/so << " new " << svi/sv << " delta " << (svi/sv  - soi/so )/dt << " ~ " <<  val << endl; 
   err += (abs((svi/sv  - soi/so )/dt - val)> 0.2);
}   //EOM

Verif(x,1)
Verif(y,2)
Verif(z,3)     
assert(err==0);
```

| The solution at $y=0.5$ |
|-------------------------|
|![][_cut]                |

[_cut]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/convect-3d/cut.png