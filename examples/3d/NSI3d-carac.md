---
name: NSI3d-carac
category: Fluid Mechanics
layout: 3d
---

## Time dependent imcompressible Navier-Stokes Equation in 3D solved by using characteristics

The Navier-Stokes equations are solved in a cube
The time independent Navier-Stokes equations for an incompressible fluid are
$$
\frac{\partial u}{\partial t}-\nu\Delta u +u\cdot\nabla u -\nabla p =0,~~~\nabla\cdot u=0~~~in ~ \Omega.
$$

The driven cavity problem has the upper side of the cube moving at velocity (1,0,0) and the other walls are at rest. By using
$$
\frac{\partial u}{\partial t} + +u\cdot\nabla u  \approx \frac1{\delta t}(u^{n+1}-u^n(x-\delta t u^n)
$$
the problem being linear at each time step, it is solved like a linear system.

The velocities are approximated with  $P^2$ elements and the pressure by a $P^2$ element (Taylor-Hood element)

The mesh is built, as usual by extrusion of a square vertically.
~~~freefem
newconvect=1;
real nu=0.01,dt=0.3;
real alpha=1./dt,alpha2=sqrt(alpha);

int nn=10;

mesh Th2=square(nn,nn);
fespace Vh2(Th2,P2);
Vh2 ux,uz,p2;
int[int] rup=[0,2],  rdown=[0,1], rmid=[1,1,2,1,3,1,4,1];
real zmin=0,zmax=1;

mesh3 Th=buildlayers(Th2,nn,
  zbound=[zmin,zmax],
  // region=r1, 
  labelmid=rmid, 
  reffaceup = rup,
  reffacelow = rdown);
~~~
The finite element spaces are defined and the matrix of the linear system for the Stokes system because it will be our initial start in time
~~~freefem
fespace VVh(Th,[P23d,P23d,P23d,P13d]);
fespace Vh(Th,P23d);
fespace Ph(Th,P13d);
macro Grad(u) [dx(u),dy(u),dz(u)]// EOM
macro div(u1,u2,u3) (dx(u1)+dy(u2)+dz(u3)) //EOM
  
  varf vStokes([u1,u2,u3,p],[v1,v2,v3,q]) = 
  int3d(Th,qforder=3)( Grad(u1)'*Grad(v1) +  Grad(u2)'*Grad(v2) +  Grad(u3)'*Grad(v3)
             - div(u1,u2,u3)*q - div(v1,v2,v3)*p + 1e-10*q*p ) 
 +  on(2,u1=1.,u2=0,u3= 0)
 + on(1,u1=0,u2=0,u3=0)
 ;

cout << "b  mat " << endl;

matrix A=vStokes(VVh,VVh);
cout << "e  mat " << endl;
~~~
The solution method is defined and the right hand side $b$.
~~~freefem
set(A,solver=sparsesolver,dimKrylov=1000);
cout << "e fac  mat " << endl;
real[int] b= vStokes(0,VVh);

VVh [u1,u2,u3,p];
VVh [X1,X2,X3,Xp];
VVh [x1,x2,x3,xp]=[x,y,z,0];

u1[]= A^-1 * b;  // Solves the Stokes problem

ux= u1(x,0.5,y);
uz= u3(x,0.5,y);
p2= p(x,0.5,y);
plot([ux,uz],p2,cmm=" cut y = 0.5",wait=1);
~~~
Now the matrix for the Navier-Stokes equations is defined
~~~freefem
macro XX1() (x-u1*dt)//
macro XX2() (y-u2*dt)//
macro XX3() (z-u3*dt)//

  varf vNS([uu1,uu2,uu3,p],[v1,v2,v3,q]) = 
  int3d(Th)( alpha*(uu1*v1+uu2*v2+uu3*v3) + nu*(Grad(uu1)'*Grad(v1) +  Grad(uu2)'*Grad(v2) +  Grad(uu3)'*Grad(v3)) //'
  - div(uu1,uu2,uu3)*q - div(v1,v2,v3)*p + 1e-10*q*p ) 
  + on(2,uu1=1,uu2=0,uu3=0)
  + on(1,uu1=0,uu2=0,uu3=0)
 
    +  int3d(Th,optimize=1,qforder=4)(   alpha*(  convect([u1,u2,u3],-dt,u1)*v1  +   convect([u1,u2,u3],-dt,u2)*v2  +   convect([u1,u2,u3],-dt,u3)*v3 )  ) ;

cout << " build  A" << endl;
A = vNS(VVh,VVh);
cout << " fac A" << endl;

set(A,solver=sparsesolver,dimKrylov=1000);
~~~
The time loop is here
~~~freefem
real t=0;
for(int i=0;i<10;++i)
  {
    t += dt;
    cout << " iteration " << i << " t = " << t << endl;
    X1[]=x1[]+u1[]*(-dt);
    b=vNS(0,VVh);
    verbosity=2;
    u1[]= A^-1 * b;
    ux= u1(x,0.5,y);
    uz= u3(x,0.5,y);
    p2= p(x,0.5,y);
    plot([ux,uz],p2,cmm=" cut y = 0.5, time ="+t,wait=0);
    // to store the solution one uses these
    if(i%5==6)
    {
      exec("mkdir dd");
      string prefu="dd/u-"+(100+i);
      string prefp="dd/p-"+(100+i);
      savemesh(Th,prefu+".mesh");
      savemesh(Th,prefp+".mesh");
     
      ofstream file(prefu+".bb"); 
      ofstream filep(prefp+".bb"); 
      Ph up1=u1,up2=u2,up3=u3,pp=p;
      file << "3 1 3 "<< up1[].n << " 2 \n";
      filep << "3 1 1 "<< pp[].n << " 2 \n";
      for (int j=0;j<up1[].n ; j++)  
	{
	  file << up1[][j] <<" " <<up2[][j] <<" "<< up3[][j] <<"\n";
	  filep << pp[][j] <<  endl; 
	}  
    }
  }
plot([ux,uz],p2,cmm=" cut y = 0.5, time ="+t,wait=1);
~~~

|The solution with a  1000 vertices mesh |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/NSI3d-carac/solution.png