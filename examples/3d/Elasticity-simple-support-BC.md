---
name: Elasticity-simple-support-BC
category: Elasticity
folder: 3d
---
## Linear elasticity with clamped parts of boundary

Homogeneous Dirichlet Boundary conditions are imposed on 2 Edge in 3d for the linear elasticity PDEs.
 The two edges are at the intersection of faces with label 2,5, x==1 (resp.  4,5, x==0).
 The deformation of the structure is due to its own weight.
~~~freefem
mesh3 Th=cube(20,5,5,[x*4,y,z]);
fespace Wh(Th,[P2,P2,P2]) ;

real E = 21.5e4;
real sigma = 0.29;
real mu = E/(2*(1+sigma));
real lambda = E*sigma/((1+sigma)*(1-2*sigma));
real gravity = -100;
~~~
The equations of linear elasticity are used
~~~freefem
real sqrt2=sqrt(2.);
macro epsilon(u1,u2,u3)  [dx(u1),dy(u2),dz(u3),(dz(u2)+dy(u3))/sqrt2,(dz(u1)+dx(u3))/sqrt2,(dy(u1)+dx(u2))/sqrt2] // EOM
macro div(u1,u2,u3) ( dx(u1)+dy(u2)+dz(u3) ) // EOM
Wh [u1,u2,u3],[v1,v2,v3];
varf von5([u1,u2,u3],[v1,v2,v3]) = on(5,u1=1,u2=1,u3=1);
varf von24([u1,u2,u3],[v1,v2,v3]) = on(2,4,u1=1,u2=1,u3=1);

Wh [au1,au2,au3];
{
	real[int] w5=von5(0,Wh, tgv=1); //  find dof on face 5 
	real[int] w24=von24(0,Wh, tgv=1); //  find dof on face 2 and 4 
	au1[] = w5 .* w24; // 1 do intersect  face 5 and 2 4 
	plot(au3,wait=1); // see for debugging 
}
// so array au1[] is non zero on dof on edges  (2,5) and (4,5)

varf  Lame([u1,u2,u3],[v1,v2,v3]) =
  int3d(Th)(  
	     lambda*div(u1,u2,u3)*div(v1,v2,v3)	
	    +2.*mu*( epsilon(u1,u2,u3)'*epsilon(v1,v2,v3) ) //')
	      )
  + int3d(Th) (gravity*v3)
;
~~~
The solution is computed
~~~freefem
matrix A = Lame(Wh,Wh,sym=1,positive=1,solver="CG");
cout << " half "<< A.half << " nnz "<< A.nnz << endl;
  real[int] b= Lame(0,Wh); 
  //  put tgv = -2 on matrix A 
  setBC(A,au1[],-2);   // 1 on diagonal and  0 on row  i et column i if  au1[][i] !=0  
  b =  au1[] ? 0 : b;
  set(A,solver= "CHOLMOD"); 
  cout << A.nnz << endl; 
   u1[] = A^-1*b;
~~~
The rest is for visualization, in particular the initial mesh is moved with the displacements just computed
~~~freefem
cout << " || u|| =" << u1[].linfty << endl;
 real cc = 0.5/u1[].linfty;
 mesh3 Thm = movemesh(Th,[x+u1*cc,y+u2*cc,z+u3*cc]) ;
 plot(Th,Thm,wait=1);
~~~

| The initial mesh and the deformed mesh |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/Elasticity-simple-support-BC/solution.png