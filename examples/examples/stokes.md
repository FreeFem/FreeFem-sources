---
name: Stokes
category: fluid
layout: example
---

# The Stokes system for creeping flow
The system is set in 2 dimensions:
$$
-\Delta u + \nabla p =\vec f, \quad
 \nabla\cdot u = 0   \hbox{ in }\Omega ,~~
 u _{|\Gamma}= \vec g
 $$
with $\Omega$ the  unit square, $\vec f=[0,0]^T$ and $\vec g={\bf 1}_{\Gamma_3}[1,0]^T$.

A possible variational formulation is
$$∀(v,q),~~∫_Ω ∇u:∇v−∫_Ωp\nabla\cdot v=0,
−∫_Ω\nabla\cdot u q−∫_Ωϵpq=0.
$$
The velocity is approximated with the $P^1+bubble$ element and the pressure by $P^1$ continuous functions.
A small penalization parameter is added to the formuation to secure uniqueness.
~~~freefem
int n = 3; // mesh quality
mesh Th = square(10*n, 10*n);

// Fespace
fespace Uh(Th, P1b);
Uh u, v;
Uh uu, vv;
fespace Ph(Th, P1);
Ph p, pp;

// Problem
solve stokes([u, v, p], [uu, vv, pp])
  = int2d(Th)(
      dx(u)*dx(uu) + dy(u)*dy(uu)
    + dx(v)*dx(vv) + dy(v)*dy(vv)
    + dx(p)*uu + dy(p)*vv
    + pp*(dx(u) + dy(v))
    -1e-10*p*pp
  )
  + on(1, 2, 4, u=0, v=0)
  + on(3, u=1, v=0)
  ;

// Plot
plot([u,v],p,wait=1);
~~~
## Results

| Isovalue lines of the pressure and vector display of the velocity |
| --------------                                                    |
| ![][_solution]                                                    |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/examples/stokes/solution.png
