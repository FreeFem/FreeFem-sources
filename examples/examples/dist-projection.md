---
name: dist-projection
category: Tools
layout: example
---

# Test the function "dist"
## and the projection on a curve built with meshL
The function $\texttt{dis(C)}$ returns the distance from the current point to $C$. Let us build a triangulation of a curve (here a cercle)  $C$, i.e. cut the circle into 20 segments.
~~~freefem
load "msh3"
border C(t=0,2*pi){ x= cos(t); y=sin(t); label=1; region =1;}
meshL Lh =buildmeshL(C(20));
~~~
Let us build a finite element space of $P^1$ functions on a square centered at (0,0) of size 3$\times$3.
~~~freefem
mesh Th= square(20,20,[(x-0.5)*3,(y-0.5)*3]);
fespace Vh(Th,P1);
Vh d=dist(Lh);
~~~
The last line defines a $P^1$ function on the square which is equal to the distance to $C$ of any point in the square Th.
Next, $U=[U_x,U_y]^T$ is the vector which starts from $P=[x,y]^T$ and ends at the projection of $P$ on $C$.
~~~freefem
Vh Ux = projection(Lh).x-x;
Vh Uy = projection(Lh).y-y;
plot(d, [Ux,Uy],Lh, wait=1);
~~~
The plot displays the level curves of $d$ in 3D and $U$.
![][_Ud]
Now the same is done with another triangulation of a square, Th2. Parameters nu containes the triangle number of the projection and ph the projection point. So Th2[nu][0] is the first vertex of triangle nu.
~~~freefem
int nu;
R3 ph;
mesh Th2= square(10,10,[(x-0.5)*3,(y-1.06)*3]);
Ux = projection(Th2,nu=nu,Phat=ph).x-x;
Uy = projection(Th2,nu=nu,Phat=ph).y-y;
x=0;y=0;z=0;
 cout << projection(Th2,nu=nu,Phat=ph) <<", " << nu <<","<< ph << endl; 
 cout << Th2[nu][0].x << endl;
 cout << Th2[nu][0].y  << endl;
 cout << Th2[nu][1].x << endl;
 cout << Th2[nu][1].y  << endl;
 cout << Th2[nu][2].x << endl;
 cout << Th2[nu][2].y  << endl;
plot([Ux,Uy],Th2, wait=1);
~~~

| The function d to Th2 |
| --------------------- |
| ![][_plot2]           |

Finally some more advanced operations are done in 3D: signeddist(ThS) is the distance from points in the cube to the square ThS.
~~~freefem
load "msh3";
//meshL Lh2 = extract(Th2);
meshS ThS= square3(10,10,[(x-0.5)*3,(y-1.06)*3,(x+y)/2]);
mesh3 Th3=cube(10,10,10,[(x-0.5)*3,(y-1.06)*3,z*2-1]);
fespace Uh(Th3,P1);
Uh d3=signeddist(ThS);
plot(d3,wait=1);
~~~

| The signeddist function |
| ----------------------- |
| ![][_plot3]             |

[_Ud]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/examples/dist-projection/Ud.png

[_plot2]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/examples/dist-projection/plot2.png

[_plot3]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/examples/dist-projection/plot3.png
