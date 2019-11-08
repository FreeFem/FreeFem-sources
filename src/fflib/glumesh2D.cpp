// $Id$

#include <iostream>
#include <cfloat>
#include <cmath>
#include <complex>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"


#include "FESpacen.hpp"
#include "FESpace.hpp"
#include "HashMatrix.hpp"

#include "SparseLinearSolver.hpp"

#include "MeshPoint.hpp"
#include "Operator.hpp"
#include "lex.hpp"

#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"

#include <set>
#include <vector>
#include <fstream>


using namespace  Fem2D;

class listMesh {
public:
  list<Mesh const  *> *lth;
  void init()  { lth=new list<Mesh const  *>;}
  void destroy() { delete lth;}
  listMesh(Stack s,Mesh const*th) : lth(Add2StackOfPtr2Free(s,new list<Mesh const*>)) { lth->push_back(th);}
  listMesh(Stack s,Mesh const*const tha,Mesh const* const thb) : lth(Add2StackOfPtr2Free(s,new list<Mesh const*>)) { lth->push_back(tha);lth->push_back(thb);}
  listMesh(Stack s,const listMesh &l,Mesh const*const th) : lth(Add2StackOfPtr2Free(s,new list<Mesh const*>(*l.lth))) { lth->push_back(th);}

};


Mesh * GluMesh(list<Mesh const *> const & lth, long labtodel = -1)
{
  int nbt=0;
  int neb=0;
  int nebx=0;
  int nbv=0;
  int nbvx=0;
  double hmin=1e100;
  R2 Pn(1e100,1e100),Px(-1e100,-1e100);
  const  Mesh * th0=0;
    int kk=0;
  for(list<Mesh const *>::const_iterator i=lth.begin();i != lth.end();++i)
    {
       if(! *i ) continue;
        ++kk;
       const Mesh &Th(**i);
      th0=&Th;
      if(verbosity>1)  cout << " GluMesh + "<< Th.nv << " " << Th.nt << endl;
      nbt+= Th.nt;
      nbvx += Th.nv;
      nebx += Th.neb;
      for (int k=0;k<Th.nt;k++)
	for (int e=0;e<3;e++)
	  hmin=min(hmin,Th[k].lenEdge(e));

      for (int i=0;i<Th.nv;i++)
	{
	  R2 P(Th(i));
	  Pn=Minc(P,Pn);
	  Px=Maxc(P,Px);
	}
    }
    if(kk==0) return 0; //  no mesh ...
  if(verbosity>2)
    cout << "      - hmin =" <<  hmin << " ,  Bounding Box: " << Pn
	 << " "<< Px << endl;

  Vertex * v= new Vertex[nbvx];
  Triangle *t= new Triangle[nbt];
  Triangle *tt=t;
  BoundaryEdge *b= new BoundaryEdge[nebx];
  BoundaryEdge *bb= b;

  ffassert(hmin>Norme2(Pn-Px)/1e9);
  double hseuil =hmin/10.;


  FQuadTree *quadtree=new Fem2D::FQuadTree(th0,Pn,Px,0);
  {
    map<pair<int,int>,int> bbe;
    for(list<Mesh const  *>::const_iterator i=lth.begin();i != lth.end();++i)
      {
          if(! *i ) continue;
	const Mesh &Th(**i);
	if(!*i) continue;
	if(verbosity>1)  cout << " GluMesh + "<< Th.nv << " " << Th.nt << endl;

	for (int ii=0;ii<Th.nv;ii++)
	{
	  const Vertex &vi(Th(ii));
	  Vertex * pvi=quadtree->ToClose(vi,hseuil);
	  if(!pvi) {
	    v[nbv].x = vi.x;
	    v[nbv].y = vi.y;
	    v[nbv].lab = vi.lab;
	    quadtree->Add(v[nbv++]);
	  }
	}

	for (int k=0;k<Th.nt;k++)
	  {
	    const Triangle  &K(Th[k]);
	    int i0=quadtree->ToClose(K[0],hseuil)-v; //NearestVertex(K[0])-v;
	    int i1=quadtree->ToClose(K[1],hseuil)-v; //NearestVertex(K[1])-v;
	    int i2=quadtree->ToClose(K[2],hseuil)-v; //NearestVertex(K[2])-v;
	    (*tt++).set(v,i0,i1,i2,K.lab);
	  }


	for (int k=0;k<Th.neb;k++)
	  {
	    const BoundaryEdge & be(Th.bedges[k]);
	    int i0=quadtree->ToClose(be[0],hseuil)-v;
	    int i1=quadtree->ToClose(be[1],hseuil)-v;
	    int ii0=i0,ii1=i1;
	    if(ii1<ii0) Exchange(ii0,ii1);
	    pair<int,int> i01(ii0,ii1);
	    if(be.lab != labtodel && bbe.find(i01) == bbe.end())
	    {
	      (*bb++).set(v,i0,i1,be.lab);
	      bbe[i01]=	  neb++;
	    }

	  }

      }
  }

  delete quadtree;


  if(verbosity>1)
    {
      cout << "     Nb points : "<< nbv << " , nb edges : " << neb << endl;
      cout << "     Nb of glu point " << nbvx -nbv;
      cout << "     Nb of glu  Boundary edge " << nebx-neb;
    }

  {
    Mesh * m = new Mesh(nbv,nbt,neb,v,t,b);
    R2 Pn,Px;
    m->BoundingBox(Pn,Px);
    m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv);
    return m;
  }

}

template<class RR,class AA=RR,class BB=AA>
struct Op2_addmesh: public binary_function<AA,BB,RR> {
  static RR f(Stack s,const AA & a,const BB & b)
  { return RR(s, a, b );}
};

template<bool INIT,class RR,class AA=RR,class BB=AA>
struct Op2_setmesh: public binary_function<AA,BB,RR> {
  static RR f(Stack stack, const AA & a,const BB & b)
  {
    ffassert(a );
    pmesh  p=GluMesh(*(b.lth));

      if(!INIT &&  *a) (**a).destroy() ;
    return *a=p,a;
  }
};

class SetMesh_Op : public E_F0mps
{
public:
  Expression a;

  static const int n_name_param =2+2+2+2+2; //  add nbiter FH 30/01/2007 11 -> 12
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];

  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const{ return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  long  arg(int i,Stack stack, long  a ) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
  bool   arg(int i,Stack stack, bool   a ) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}


public:
  SetMesh_Op(const basicAC_F0 &  args,Expression aa) : a(aa) {
    args.SetNameParam(n_name_param,name_param,nargs);
      if( nargs[0] && nargs[2] )
	  CompileError("uncompatible change (Th, label= , refe=  ");
      if( nargs[1] && nargs[3] )
	  CompileError("uncompatible change (Th, region= , reft=  ");

  }

  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type SetMesh_Op::name_param[]= {
  {  "refe", &typeid(KN_<long> )},
  {  "reft", &typeid(KN_<long> )},
  {  "label", &typeid(KN_<long> )},
  {  "region", &typeid(KN_<long> )},
  {  "renumv",&typeid(KN_<long>)},
  {  "renumt",&typeid(KN_<long>)},
  {  "flabel", &typeid(long)},
  {  "fregion", &typeid(long)},
  {  "rmledges", &typeid(long)},
  {  "rmInternalEdges", &typeid(bool)}
};

int  ChangeLab(const map<int,int> & m,int lab)
{
  map<int,int>::const_iterator i=m.find(lab);
  if(i != m.end())
    lab=i->second;
  return lab;
}

AnyType SetMesh_Op::operator()(Stack stack)  const
{
    MeshPoint *mp=MeshPointStack(stack),smp=*mp;
  Mesh * pTh= GetAny<Mesh *>((*a)(stack));
  Mesh & Th=*pTh;
  Mesh *m= pTh;
  int nbv=Th.nv; // nombre de sommet
  int nbt=Th.nt; // nombre de triangles
  int neb=Th.neb; // nombre d'aretes fontiere
  KN<long> zz;
  KN<long> nre (arg(0,stack,arg(2,stack,zz)));
  KN<long> nrt (arg(1,stack,arg(3,stack,zz)));
  KN<long> rv (arg(4,stack,zz));
  KN<long> rt (arg(5,stack,zz));
  Expression flab = nargs[6] ;
  Expression freg = nargs[7] ;
  bool  rm_edge = nargs[8];
  long  rmlabedges (arg(8,stack,0L));
  bool  rm_i_edges = (arg(9,stack,false));

  bool rV =  (rv.size()== nbv);
  bool rT =  (rt.size()== nbt);
    if(verbosity>1)
	cout << "  -- SetMesh_Op: nb vertices" << nbv<< " nb Trai "<< nbt << " nb b. edges  "
	     << neb << "renum V " << rV << " , renum T "<< rT << " rm internal edges " << rm_i_edges<< endl;

  if(nre.N() <=0 && nrt.N()<=0  && !rV && ! rT  && ! flab && ! freg && !rm_i_edges &&   !rm_edge ) return m;
  ffassert( nre.N() %2 ==0);
  ffassert( nrt.N() %2 ==0);
  map<int,int> mape,mapt;
  int z00 = false;
  for(int i=0;i<nre.N();i+=2)
    { z00 = z00 || ( nre[i]==0 && nre[i+1]==0);
      if(nre[i] != nre[i+1])
	mape[nre[i]]=nre[i+1];
    }
  for(int i=0;i<nrt.N();i+=2)
    mapt[nrt[i]]=nrt[i+1];
  int nebn =0;
    for(int i=0;i<neb;++i)
      {
	int l0,l1=ChangeLab(mape,l0=m->bedges[i].lab) ;
	 nebn++;
      }

  Vertex * v= new Vertex[nbv];
  Triangle *t= new Triangle[nbt];
  BoundaryEdge *b= new BoundaryEdge[nebn];
  // generation des nouveaux sommets
  Vertex *vv=v;
  // copie des anciens sommets (remarque il n'y a pas operateur de copy des sommets)
  for (int i=0;i<nbv;i++)
    {
      int ii=rV? rv(i): i;
      vv = v + ii;
     Vertex & V=Th(i);
     vv->x=V.x;
     vv->y=V.y;
     vv->lab = V.lab;

   }

  //  generation des triangles
  int nberr=0;
  R2 PtHat(1./3.,1./3.);
  for (int i=0;i<nbt;i++)
    {
      int ii= rT ? rt(i) : i;
      int i0=Th(i,0), i1=Th(i,1),i2=Th(i,2);
      if(rV) {
	  i0=rv(i0);
	  i1=rv(i1);
	  i2=rv(i2);
      }
      // les 3 triangles par triangles origines
      t[ii].set(v,i0,i1,i2,ChangeLab(mapt,Th[i].lab));
      if(freg)
	{
	  mp->set(Th,Th[i](PtHat),PtHat,Th[i],Th[i].lab);
	  t[ii].lab =GetAny<long>( (* freg)(stack)) ;
	}

    }
  int  nrmedge=0;
  // les arete frontieres qui n'ont pas change
  BoundaryEdge * bb=b;
  for (int i=0;i<neb;i++)
    {
      int ke,k =Th.BoundaryElement(i,ke);
      int kke,kk= Th.ElementAdj(k,kke=ke);
      int intern = ! (( kk == k )|| ( kk < 0)) ;
      const   Triangle &K(Th[k]);
      int i1=Th(Th.bedges[i][0]);
      int i2=Th(Th.bedges[i][1]);

	if(rV) {
	    i1=rv(i1);
	    i2=rv(i2);
	}

      int l0,l1=ChangeLab(mape,l0=m->bedges[i].lab) ;
      mp->set(Th,Th[k](PtHat),PtHat,Th[k],l1);
      if(flab)
      {
          R2 E=K.Edge(ke);
          double le = sqrt((E,E));
          double sa=0.5,sb=1-sa;
          R2 PA(TriangleHat[VerticesOfTriangularEdge[ke][0]]),
          PB(TriangleHat[VerticesOfTriangularEdge[ke][1]]);
          R2 Pt(PA*sa+PB*sb );
          MeshPointStack(stack)->set(Th,K(Pt),Pt,K,l1,R2(E.y,-E.x)/le,ke);
	  l1 =GetAny<long>( (*flab)(stack)) ;
      }
     if( !( intern && ( rm_i_edges || (rmlabedges && (l1 == rmlabedges)) )   ))
      *bb++ = BoundaryEdge(v,i1,i2,l1);
     else ++nrmedge;
    }
    nebn -= nrmedge;
    if(nrmedge && verbosity > 2) cout << "   change  mesh2 : number of removed  internal edges " << nrmedge << endl;
  assert(nebn==bb-b);
  m =  new Mesh(nbv,nbt,nebn,v,t,b);

  R2 Pn,Px;
  m->BoundingBox(Pn,Px);
  m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv);
  Add2StackOfPtr2FreeRC(stack,m);
    *mp=smp;
  return m;
}

class SetMesh : public OneOperator { public:
typedef Mesh const *pmesh;
    SetMesh() : OneOperator(atype<pmesh>(),atype<pmesh>() ) {}

  E_F0 * code(const basicAC_F0 & args) const
  {
    return  new SetMesh_Op(args,t[0]->CastTo(args[0]));
  }
};

Mesh* GluMeshtab (KN<pmesh> *const &tab, long const &lab_delete) {
  list<Mesh const *> l;
  for (int i=0; i<tab->n; i++)
    l.push_back((*tab)[i]);
  return GluMesh(l,lab_delete);
}

struct Op_GluMeshtab: public OneOperator {
    typedef const Mesh *pmesh;
    class Op: public E_F0mps   {
    public:
        static basicAC_F0::name_and_type name_param [];
        static const int n_name_param = 1;
        Expression nargs[n_name_param];
        Expression getmeshtab;
        long arg (int i, Stack stack, long a) const {return nargs[i] ? GetAny<long>((*nargs[i])(stack)) : a;}

        Op (const basicAC_F0 &args, Expression t): getmeshtab(t)
        {args.SetNameParam(n_name_param, name_param, nargs);}

        AnyType operator () (Stack s)  const;
    };

    E_F0*code (const basicAC_F0 &args) const
    {return new Op(args, t[0]->CastTo(args[0]));}

    Op_GluMeshtab ():
    OneOperator(atype<const pmesh>(), atype<KN<pmesh> *>()) {};
};
basicAC_F0::name_and_type Op_GluMeshtab::Op::name_param[Op_GluMeshtab::Op::n_name_param] =
{
    {"labtodel", &typeid(long)}
};
AnyType Op_GluMeshtab::Op::operator () (Stack stack)  const {
    KN<const Mesh *> *tab = GetAny<KN<const Mesh *> *>((*getmeshtab)(stack));
    long labtodel = arg(0, stack, LONG_MIN);
    Mesh *Tht = GluMeshtab(tab, labtodel);

    Add2StackOfPtr2FreeRC(stack, Tht);
    return Tht;
}

//  truc pour que la fonction
// Init::Init() soit appele a moment du chargement dynamique
// du fichier
//
#ifndef  DYNAMICS_LIBS
void init_glumesh2D()
{
  Dcl_Type<listMesh>();
  typedef Mesh const  *pmesh;

  if(verbosity>2)
    cout << " glumesh2D " ;
  TheOperators->Add("+",new OneBinaryOperator_st< Op2_addmesh<listMesh,pmesh,pmesh>  >      );
  TheOperators->Add("+",new OneBinaryOperator_st< Op2_addmesh<listMesh,listMesh,pmesh>  >      );
  TheOperators->Add("=",new OneBinaryOperator_st< Op2_setmesh<false,pmesh*,pmesh*,listMesh>  >     );
  TheOperators->Add("<-",new OneBinaryOperator_st< Op2_setmesh<true,pmesh*,pmesh*,listMesh>  >     );

  Global.Add("change","(",new SetMesh);

  Global.Add("gluemesh", "(", new Op_GluMeshtab);
}
#else
class Init { public:
  Init();
};

static Init init;  //  une variable globale qui serat construite  au chargement dynamique

Init::Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  Dcl_Type<listMesh>();
  typedef Mesh *pmesh;

  if (verbosity)
    cout << "  glumesh2D " ;
  TheOperators->Add("+",new OneBinaryOperator_st< Op2_addmesh<listMesh,pmesh,pmesh>  >      );
  TheOperators->Add("+",new OneBinaryOperator_st< Op2_addmesh<listMesh,listMesh,pmesh>  >      );
  TheOperators->Add("=",new OneBinaryOperator_st< Op2_setmesh<false,pmesh*,pmesh*,listMesh>  >     );
  TheOperators->Add("<-",new OneBinaryOperator_st< Op2_setmesh<true,pmesh*,pmesh*,listMesh>  >     );

  Global.Add("change","(",new SetMesh);

  Global.Add("gluemesh", "(", new Op_GluMeshtab);
}
#endif
