// $Id$

#include  <iostream>
#include  <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;  
#include "rgraph.hpp"
#include "RNM.hpp"
#include "ff++.hpp"
#include "AFunction_ext.hpp" // [[file:../src/fflib/AFunction_ext.hpp]]
#include "msh3.hpp"
#include <fem.hpp>
#include <cmath>

  using namespace  Fem2D;

Mesh3 const * SplitMesh4(Stack stack,Fem2D::Mesh3 const * const & pTh)
{
  assert(pTh);
  const Mesh3 & Th(*pTh);  // le maillage d'origne a decoupe
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int nbe=Th.nbe; // nombre d'aretes fontiere
  // allocation des nouveaux items du maillage  
  Vertex3 * v= new Vertex3[nbv+nbt];
  Tet *t= new Tet[nbt*4];
  Triangle3 *b  = new Triangle3[nbe];
  // generation des nouveaus sommets 
  Vertex3 *vv=v;
  // copie des anciens sommets (remarque il n'y a pas operateur de copy des sommets)
  for (int i=0;i<nbv;i++)
   {
     const Vertex3 & V=Th(i);
     vv->x=V.x;
     vv->y=V.y;
     vv->z=V.z;
     vv->lab = V.lab;
     vv++;      
   }
  // generation des points barycentre de trianngles 
  for (int k=0;k<nbt;k++)
    {
      const Tet & K=Th[k];
      R3 G= ( (R3) K[0] + K[1] + K[2] + K[3] )  / 4.;
      vv->x=G.x;
      vv->y=G.y;
      vv->z=G.z;
      vv->lab = 0;
      vv++;
    }   

  //  generation des triangles 
  Tet *tt= t; 
   
  for (int i=0;i<nbt;i++)
    {
      int i0=Th(i,0), i1=Th(i,1),i2=Th(i,2),i3=Th(i,3);
      int ii = nbv + i; // numero du 
      // les 3 triangles par triangles origines 
      int ivt[4] = {ii,i1,i2,i3};
      (*tt++).set(v,ivt,Th[i].lab);
      ivt[0] = i0;
      ivt[1] = ii;
      (*tt++).set(v,ivt,Th[i].lab);
      ivt[1] = i1;
      ivt[2] = ii;
      (*tt++).set(v,ivt,Th[i].lab);
      ivt[2] = i2;
      ivt[3] = ii;
      (*tt++).set(v,ivt,Th[i].lab);
    }

  // les arete frontieres qui n'ont pas change
  Triangle3 * bb=b;
  for (int i=0;i<nbe;i++)
  {        
      const Triangle3 &K(Th.be(i));
      int ivv[3];

      ivv[0] = Th.operator()(K[0]);
      ivv[1] = Th.operator()(K[1]);
      ivv[2] = Th.operator()(K[2]);
      (bb++)->set( v, ivv, K.lab);
  }
  //  generation de la class Mesh a partir des 3 tableaux : v,t,b
  {
    Mesh3 * m = new Mesh3(nbv+nbt,nbt*4,nbe,v,t,b);
    m->BuildGTree();
   // m->decrement();
     Add2StackOfPtr2FreeRC(stack,m);
    return m;
  }
}

//  truc pour que la fonction 
// static void Load_Init() soit appele a moment du chargement dynamique
// du fichier 
//  
/*  class Init { public:
  Init();
};

$1 */

static void Load_Init(){  // le constructeur qui ajoute la fonction "splitmesh4"  a freefem++ 
    if (verbosity > 1)
        cout << " load: Split4  " << endl;
    Global.Add("splitmesh4","(",new OneOperator1s_<Mesh3 const *,Mesh3 const *>(SplitMesh4));
}
LOADFUNC(Load_Init)
