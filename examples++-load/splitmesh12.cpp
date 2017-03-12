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

Mesh3 const * SplitMesh12(Stack stack,Fem2D::Mesh3 const * const & pTh)
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
  int nbnv = 0;
  for(int k=0; k <nbt; ++k)
      for (int e = 0; e<4;++e) {
          int ee=e;
          int kk = Th.ElementAdj(k,ee);
          if( kk>k)
              nbnv++;
      }
  Vertex3 * v= new Vertex3[nbv+nbt+nbe+nbnv];
  Tet *t= new Tet[nbt*12];
  Triangle3 *b  = new Triangle3[nbe*3];
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
  for (int i=0;i<nbe;i++)
    {        
        const Triangle3 &K(Th.be(i));
        R3 G= ( (R3) K[0] + K[1] + K[2])  / 3.;
        vv->x=G.x;
        vv->y=G.y;
        vv->z=G.z;
        vv->lab = Th.be(i).lab;
        vv++;
    }

  //  generation des triangles 
  Tet *tt= t; 
   
  for (int i=0;i<nbe;i++)
    {        
        int ki;
        int k=Th.BoundaryElement(i,ki);
        const Triangle3 &K(Th.be(i));
        int i0=Th.operator()(K[0]), i1=Th.operator()(K[1]),i2=Th.operator()(K[2]);
        int ii = nbv + nbt + i; // numero du 
        int jj = nbv + k; // numero du 
        int ivt[4] = {jj,i1,i2,ii};
        (*tt++).set(v,ivt,K.lab);
        ivt[0] = i0;
        ivt[1] = jj;
        (*tt++).set(v,ivt,K.lab);
        ivt[1] = i1;
        ivt[2] = jj;
        (*tt++).set(v,ivt,K.lab);
    }
  for (int k=0;k<nbt;k++)
  {
      for (int e = 0; e<4;++e) {
          int ee=e;
          int kk = Th.ElementAdj(k,ee);
          if( kk>k){
              const Tet & K=Th[k];
              const Tet & KAdj=Th[kk];
              int ABC[3];
              {
                  int KNo = -1;
                  {
                      for(int i = 0; i < 3 && KNo == -1; ++i) {
                          int diff = 0;
                          for(int j = 0 ; j < 4; ++j) {
                              R3 d = KAdj[j] - K[i];
                              if(d.norme() > 1e-8)
                                  ++diff;
                          }
                          if(diff == 4)
                              KNo = i;
                      }
                  }
                  if(KNo == 0)
                      ABC[0] = 1, ABC[1] = 2, ABC[2] = 3;
                  else if(KNo == 1)
                      ABC[0] = 0, ABC[1] = 2, ABC[2] = 3;
                  else if(KNo == 2)
                      ABC[0] = 0, ABC[1] = 1, ABC[2] = 3;
                  else
                      ABC[0] = 0, ABC[1] = 1, ABC[2] = 2;
                  R3 uuu = K[ABC[1]] - K[ABC[0]], vvv = K[ABC[2]] - K[ABC[0]];
                  R3 n = uuu^vvv;
                  Vertex3 dir(v[nbv+k] - v[nbv+kk]);
                  Vertex3 w0(v[nbv+kk] - K[ABC[0]]);
                  double aa = -(n,w0);
                  double bb =  (n,dir);
                  double rr =  aa / bb;
                  vv->x=v[nbv+kk].x + rr * dir.x;
                  vv->y=v[nbv+kk].y + rr * dir.y;
                  vv->z=v[nbv+kk].z + rr * dir.z;
                  vv->lab = 0;
              }
              int i0=Th.operator()(K[ABC[0]]), i1=Th.operator()(K[ABC[1]]),i2=Th.operator()(K[ABC[2]]);
              for(int ij = 0; ij < 2; ++ij) {
                  int ivt[4] = {ij==0?nbv+k:nbv+kk,i1,i2,vv - v};
                  if(det(v[ivt[0]], v[ivt[1]], v[ivt[3]]) > 0)
                      std::swap(ivt[0], ivt[3]);
                  int lab = ij==0?K.lab:KAdj.lab;
                  (*tt++).set(v,ivt,lab);
                  ivt[1] = ivt[0];
                  ivt[0] = i0;
                  (*tt++).set(v,ivt,lab);
                  ivt[2] = ivt[1];
                  ivt[1] = i1;
                  (*tt++).set(v,ivt,lab);
              }
              vv++;
          }
      }
  }

  // les arete frontieres qui n'ont pas change
  Triangle3 * bb=b;
  for (int i=0;i<nbe;i++)
  {        
      const Triangle3 &K(Th.be(i));
      int orig[3];
      int ivv[3];

      orig[0] = Th.operator()(K[0]);
      orig[1] = Th.operator()(K[1]);
      orig[2] = Th.operator()(K[2]);
      ivv[0] = nbv+nbt+i;
      ivv[1] = orig[1];
      ivv[2] = orig[2];
      (bb++)->set( v, ivv, K.lab);
      ivv[1] = ivv[0];
      ivv[0] = orig[0];
      (bb++)->set( v, ivv, K.lab);
      ivv[2] = ivv[1];
      ivv[1] = orig[1];
      (bb++)->set( v, ivv, K.lab);
  }
  //  generation de la class Mesh a partir des 3 tableaux : v,t,b
  {
    Mesh3 * m = new Mesh3(nbv+nbt+nbe+nbnv,nbt*12,nbe*3,v,t,b);
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

static void Load_Init(){  // le constructeur qui ajoute la fonction "splitmesh12"  a freefem++ 
    if (verbosity > 1)
        cout << " load: Split12  " << endl;
    Global.Add("splitmesh12","(",new OneOperator1s_<Mesh3 const *,Mesh3 const *>(SplitMesh12));
}
LOADFUNC(Load_Init)
