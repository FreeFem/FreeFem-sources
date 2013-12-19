*** /home/alh/ff/draft/src/femlib/fem.hpp	2013-12-19 11:32:21.039226795 +0100
--- /tmp/fem.hpp~other.9RdaPG	2013-12-19 11:32:21.035226766 +0100
***************
*** 24,45 ****
  
  //#include "ufunction.hpp" 
  #include "ufunction.hpp" 
! inline double norm(double x){return x*x;} 
! inline float norm(float x){return x*x;}
  #include <utility>
  #include <algorithm>
! 
! // ALH - 21/10/13 - R1.hpp needs cmath and it does include it, but since R1.hpp is only included as part of namespace
! // Fem2D, we need to make sure that cmath is called as part of the default namespace, otherwise we will get "error:
! // ‘::acos’ has not been declared" and such.
! 
! #include <cmath>
  
  // definition R
  namespace Fem2D 
  {
! 
!   // ALH - These include files are located inside the namespace definition on purpose?
  
  #include "R1.hpp"
  #include "R2.hpp"
--- 24,40 ----
  
  //#include "ufunction.hpp" 
  #include "ufunction.hpp" 
! 
  #include <utility>
  #include <algorithm>
! #include <complex>
  
  // definition R
  namespace Fem2D 
  {
! inline double norm(double x){return x*x;} 
! inline float norm(float x){return x*x;}
! template<class T> T  norm(const complex<T> &x){return std::norm(x);}
  
  #include "R1.hpp"
  #include "R2.hpp"
***************
*** 565,574 ****
    int *BoundaryAdjacencesHead;
    int *BoundaryAdjacencesLink; 
    int *TriangleConteningVertex;       
! 
!   // <<no_mesh_copy>> the copy constructor for Mesh is kept private on purpose
    Mesh(const Mesh &);
- 
    void operator=(const Mesh &);       
  };
  
--- 560,567 ----
    int *BoundaryAdjacencesHead;
    int *BoundaryAdjacencesLink; 
    int *TriangleConteningVertex;       
!   // no copy
    Mesh(const Mesh &);
    void operator=(const Mesh &);       
  };
  
***************
*** 640,646 ****
  }
  
   inline   int numSubTVertex(int N,int i,int j)
!     {  //  i,j  coordonne barycentre * N dans l'eleùent de reference.
  	i=i+j; // numerotation / diag  
  	// i,j 
  	assert(j<=i && 0<= j); 
--- 633,639 ----
  }
  
   inline   int numSubTVertex(int N,int i,int j)
!     {  //  i,j  coordonne barycentre * N dans l'ele�ent de reference.
  	i=i+j; // numerotation / diag  
  	// i,j 
  	assert(j<=i && 0<= j); 
