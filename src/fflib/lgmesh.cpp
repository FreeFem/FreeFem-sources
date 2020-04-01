// -*- Mode : c++ -*-
//
// SUMMARY  :
// USAGE    :
// ORG      :
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

/*

 This file is part of Freefem++

 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.


 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "ff++.hpp"
#include "AFunction_ext.hpp"


#include "lgmesh.hpp"

using Fem2D::Mesh;
using Fem2D::MeshPoint;

extern bool NoWait;

typedef Mesh const * pmesh;

class classBuildMesh :  public E_F0mps { public:

   typedef pmesh  Result;

   static basicAC_F0::name_and_type name_param[] ;
   static const int n_name_param =5;

    Expression nargs[n_name_param];

   Expression getborders;

   long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
   bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
   double arg(int i,Stack stack,double a) const{ return nargs[i] ? GetAny<double>( (*nargs[i])(stack) ): a;}
   KNM<double>* arg(int i,Stack stack,KNM<double>* p) const{ return nargs[i] ? GetAny<KNM<double>*>( (*nargs[i])(stack) ): p;}

    classBuildMesh(const basicAC_F0 & args)
    {
      args.SetNameParam(n_name_param,name_param,nargs);
      getborders=to<const E_BorderN *>(args[0]);
     }

    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<const E_BorderN *>());}
    static  E_F0 * f(const basicAC_F0 & args){ return new classBuildMesh(args);}
    AnyType operator()(Stack s) const ;
    operator aType () const { return atype<Result>();}

};

class classBuildMeshArray :  public E_F0mps { public:

    typedef pmesh  Result;

    //static basicAC_F0::name_and_type name_param[] ;
    static const int n_name_param =0;

//Expression nargs[n_name_param];

    Expression env,ent,enbe;
    Expression ev,et,ebe;




     classBuildMeshArray(const basicAC_F0 & args)
    {
	args.SetNameParam(0,0,0);
	env=to<long>(args[0]);
	ent=to<long>(args[1]);
	enbe=to<long>(args[2]);
	ev=to<KNM_<double> >(args[3]);
	et=to<KNM_<long> >(args[4]);
	ebe=to<KNM_<long> >(args[5]);

    }

    static ArrayOfaType  typeargs() {

	aType ffint =atype< long >();
	aType ffMi =atype< KNM_<long> >();
	aType ffMr =atype< KNM_<double> >();
	return  ArrayOfaType(ffint,ffint,ffint,ffMr,ffMi,ffMi);}

    static  E_F0 * f(const basicAC_F0 & args){ return new classBuildMeshArray(args);}
    AnyType operator()(Stack s) const {
	long nv = GetAny<long >((*env)(s));
	long nt = GetAny<long >((*ent)(s));
	long nbe = GetAny<long >((*enbe)(s));
	KNM_<double> xyl( GetAny<KNM_<double> >((*ev)(s)));
	KNM_<long> nut(	GetAny<KNM_<long> >((*et)(s)));
	KNM_<long> nube( GetAny<KNM_<long> >((*ebe)(s)));
	using  Fem2D::Vertex;
	using  Fem2D::R2;
	using  Fem2D::BoundaryEdge;
	using  Fem2D::Mesh;
	using  Fem2D::MeshPointStack;

	cout << nv << " " << nt << " " << nbe << endl;
	cout << xyl.N() << " "<< xyl.M() << endl;
	cout << nut.N() << " "<< nut.M() << endl;
	cout << nube.N() << " "<< nube.M() << endl;
	ffassert(xyl.N() >=nv && xyl.M() >= 3);
	ffassert(nut.N() >=nt && nut.M() >= 4);
	ffassert(nube.N() >=nbe && nube.M() >= 3);
	Vertex *v= new Vertex [nv];
	Triangle *t= new Triangle[nt];
	BoundaryEdge *b = new BoundaryEdge[nbe];

	for(int i=0;i<nv;++i)
	  {
	    v[i].x = xyl(i,0);
	    v[i].y = xyl(i,1);
	    v[i].lab = xyl(i,2);
	  }
	for(int i=0;i<nt;++i)
	    t[i].set(v,nut(i,0),nut(i,1),nut(i,2),nut(i,3));
	for(int i=0;i<nbe;++i)
	    b[i].set(v,nube(i,0),nube(i,1),nube(i,2));


	Mesh * pth= new Mesh(nv,nt,nbe,v,t,b);
	if(verbosity) cout << "  -- BuildMesh " << pth << " " << nv << " " << nt << " " << nbe << endl;
	R2 Pn,Px;
	pth->BoundingBox(Pn,Px);
	if(!pth->quadtree)
	    pth->quadtree=new Fem2D::FQuadTree(pth,Pn,Px,pth->nv);


	return Add2StackOfPtr2FreeRC(s,pth);//  07/2008 FH


    }
    operator aType () const { return atype<Result>();}

};


basicAC_F0::name_and_type  classBuildMesh::name_param[]= {
    {  "nbvx", &typeid(long)} ,
    {"fixeborder", &typeid(bool)},// obsolete
    {"points", &typeid(KNM<double>*)},
    {"fixedborder", &typeid(bool)},
     {"alea", &typeid(double)}
};
// modif aout 2007
class BuildMeshFile :  public E_F0mps { public:

    typedef pmesh  Result;

    static basicAC_F0::name_and_type name_param[] ;
    static const int n_name_param =1;

    Expression nargs[n_name_param];

    Expression getfilename;

    long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}

    BuildMeshFile(const basicAC_F0 & args)
    {
	args.SetNameParam(n_name_param,name_param,nargs);
	getfilename=to<string* >(args[0]);
    }

    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<string *>());}
    static  E_F0 * f(const basicAC_F0 & args){ return new BuildMeshFile(args);}
    AnyType operator()(Stack s) const ;
    operator aType () const { return atype<Result>();}

};

basicAC_F0::name_and_type  BuildMeshFile::name_param[]= {
    {  "nbvx", &typeid(long) }
};
// fin modif 2007

class MoveMesh :  public E_F0mps { public:

   typedef pmesh  Result;
   Expression getmesh;
   Expression U,V;
   int nbsol;
    vector<Expression> sol;

    MoveMesh(const basicAC_F0 & args) :nbsol(args.size()-2),sol(args.size()-2)
    {
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
 // using  Fem2D::R;
  using  Fem2D::MeshPointStack;
      args.SetNameParam();
      getmesh=to<pmesh>(args[0]);
      const E_Array * a = dynamic_cast<const E_Array *>(args[1].LeftValue());

      ffassert(a);
      if (a->size() !=2) CompileError("movemesh(Th,[u,v],...) need 2 componate in array ",atype<pmesh>());
      U=to<double>( (*a)[0]);
      V=to<double>( (*a)[1]);

      for (int i=2;i<args.size();i++)
        sol[i-2]=to<double>(args[i]);
    }
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),atype<E_Array>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new MoveMesh(args);}
    AnyType operator()(Stack s) const ;
      operator aType () const { return atype<Result>();}

};

class SplitMesh :  public E_F0mps { public:

   typedef pmesh  Result;
   Expression getmesh;
   Expression U;

    SplitMesh(const basicAC_F0 & args)
    {
      args.SetNameParam();
      getmesh=to<pmesh>(args[0]);
      U = to<long>(args[1]);
     }
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),atype<long>());}
    static  E_F0 * f(const basicAC_F0 & args){ return new SplitMesh(args);}
    AnyType operator()(Stack s) const ;
    operator aType () const { return atype<Result>();}

};

class SaveMesh :  public E_F0 { public:

   typedef pmesh  Result;
   Expression getmesh;
   Expression filename;
   Expression xx,yy,zz;
   SaveMesh(const basicAC_F0 & args)
    {
      xx=0;
      yy=0;
      zz=0;
      args.SetNameParam();
      getmesh=to<pmesh>(args[0]);
      filename=to<string*>(args[1]);
      if (args.size() >2)
        {
          const E_Array * a = dynamic_cast<const E_Array *>(args[2].LeftValue());
          if (!a) CompileError("savemesh(Th,\"filename\",[u,v,w],...");
          int k=a->size() ;
         // cout << k << endl;
          if ( k!=2 && k !=3) CompileError("savemesh(Th,\"filename\",[u,v,w]) need 2 or 3  componate in array ",atype<pmesh>());
          xx=to<double>( (*a)[0]);
          yy=to<double>( (*a)[1]);
          if(k==3)
           zz=to<double>( (*a)[2]);
         }

   }
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),atype<string*>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new SaveMesh(args);}
    AnyType operator()(Stack s) const ;

};


class Adaptation :   public E_F0mps { public:
  typedef pmesh  Result;

  static basicAC_F0::name_and_type name_param[] ;
  static const int n_name_param =28;

  int nbsol;
  Expression nargs[n_name_param];
  Expression getmesh;
  Expression em11,em22,em12;
  int  typesol[100];
  vector<Expression> sol;
  int nbcperiodic;
  Expression *periodic;


  double arg(int i,Stack stack,double a) const { return nargs[i] ? GetAny<double>( (*nargs[i])(stack) ): a;}
  long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
  long* arg(int i,Stack stack,long* a) const{ return nargs[i] ? GetAny<long*>( (*nargs[i])(stack) ): a;}
  bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
  int arg(int i,Stack stack,int a) const{ return nargs[i] ? GetAny<int>( (*nargs[i])(stack) ): a;}

  Adaptation(const basicAC_F0 & args) :nbsol(args.size()-1),sol(args.size()-1)
  {
      em11=0;
      em22=0;
      em12=0;

      args.SetNameParam(n_name_param,name_param,nargs);
      getmesh=to<pmesh>(args[0]);
      int ksol=0;
      ffassert(nbsol<100);
      for (int i=1;i<nbsol+1;i++)
         if (args[i].left()==atype<E_Array>())
          {
            const E_Array * a = dynamic_cast<const E_Array *>(args[i].LeftValue());
            ffassert(a);
            ksol+=a->size();
          }
         else
           ksol++;
      sol.resize(ksol);
      ksol=0;
      for (int i=1;i<nbsol+1;i++)
         if (args[i].left()==atype<E_Array>())
          {
            const E_Array * a = dynamic_cast<const E_Array *>(args[i].LeftValue());
            ffassert(a);
             int N=a->size();
            typesol[i-1]=N-1; // ok en 2D
            if (N<=4) {
             for (int j=0;j<N;j++)
              sol[ksol++]=to<double>((*a)[j]); }
            else {
             lgerror(" Adaptation vecteur a plus de 4 compossantes inconnue");
            }

          }
         else {
          typesol[i-1]=0;
          sol[ksol++]=to<double>(args[i]);
          }
    const E_Array * expmetrix = dynamic_cast<const E_Array *>(nargs[24]);

  if(expmetrix)
   {
      if(expmetrix->nbitem()!=3)
        ExecError("\nSorry we wait an array with 3 componants in: metrix=[m11,m12,m22]");

      em11=(*expmetrix)[0];
      em12=(*expmetrix)[1];
      em22=(*expmetrix)[2];

      if(  (*expmetrix)[0].left()!= atype<KN<double> *>() )
          CompileError("Sorry the fist  array componant in metrix=[m11,m12,m22] must be vector");
      if( (*expmetrix)[1].left()!= atype<KN<double> *>() )
          CompileError("Sorry the second  array componant in metrix=[m11,m12,m22] must be vector");
      if( (*expmetrix)[2].left()!= atype<KN<double> *>() )
          CompileError("Sorry the third  array componant in metrix=[m11,m12,m22] must be vector");


     }
     nbcperiodic=0;
     periodic=0;
     GetPeriodic(2,nargs[25],nbcperiodic,periodic);

   }

    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new Adaptation(args);}
    AnyType operator()(Stack s) const ;
    operator aType () const { return atype<pmesh>();}
};


 basicAC_F0::name_and_type Adaptation::name_param[Adaptation::n_name_param] = {
       {   "hmin",             &typeid(double)},  // �
       {   "hmax",             &typeid(double)},
       {   "err",              &typeid(double)},
       {   "errg",             &typeid(double)},
       {   "nbvx",             &typeid(long)},        // 4
       {   "nbsmooth",         &typeid(long)},
       {   "nbjacoby",         &typeid(long)},
       {   "ratio",            &typeid(double)},
       {   "omega",            &typeid(double)},
       {   "iso",              &typeid(bool)},         // 9
       {   "abserror",         &typeid(bool)},
       {   "cutoff",           &typeid(double)},
       {   "verbosity",        &typeid(long)},
       {   "inquire",          &typeid(bool)},
       {   "splitpbedge",      &typeid(bool)}, // 14
       {   "maxsubdiv",        &typeid(double)},
       {   "anisomax",         &typeid(double)},
       {   "rescaling",        &typeid(bool)},
       {   "keepbackvertices", &typeid(bool)},
       {   "IsMetric",         &typeid(bool)},    // 19
       {   "power",            &typeid(double)},    // 20
       {   "thetamax",         &typeid(double)},
       {   "splitin2",         &typeid(bool) },
       {  "nomeshgeneration", &typeid(bool) },
       {   "metric"          ,  &typeid(E_Array)},  // 24
       {   "periodic"        ,  &typeid(E_Array) },// 25
       { "requirededges",    &typeid(KN_<long> ) }, // 26
       { "warning",    &typeid(long *) } // 27

    };

struct Op_trunc_mesh : public OneOperator {

    class Op: public E_F0mps   { public:
      static basicAC_F0::name_and_type name_param[] ;
      static const int n_name_param =9;
      Expression nargs[n_name_param];

      Expression getmesh,bbb;
      long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
      bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}

       KN<long> *  arg(int i,Stack stack) const{ return nargs[i] ? GetAny<KN<long> *>( (*nargs[i])(stack) ): 0;}
      Op(const basicAC_F0 &  args,Expression t,Expression b) : getmesh(t),bbb(b)
        { args.SetNameParam(n_name_param,name_param,nargs); }
      AnyType operator()(Stack s)  const ;
     };
    E_F0 * code(const basicAC_F0 & args) const
     { return new Op(args,to<pmesh>(args[0]),to<bool>(args[1])) ;}
   Op_trunc_mesh() :
     OneOperator(atype<pmesh>(),atype<pmesh>(),atype<bool>()) {}
};

basicAC_F0::name_and_type Op_trunc_mesh::Op::name_param[Op_trunc_mesh::Op::n_name_param] =
 {
   {  "split",             &typeid(long)},
   {  "label",             &typeid(long)},
   { "new2old", &typeid(KN<long>*)},  //  ajout FH pour P. Jolivet jan 2014
   { "old2new", &typeid(KN<long>*)},   //  ajout FH pour P. JoLivet jan 2014
     { "renum",&typeid(bool)},
//  add jan 2019
     {  "labels", &typeid(KN_<long> )},
     {  "region", &typeid(KN_<long> )},
     {  "flabel", &typeid(long)},
     {  "fregion", &typeid(long)},
 };

AnyType classBuildMesh::operator()(Stack stack)  const {
    const E_BorderN * borders = GetAny<const E_BorderN *>((*getborders)(stack));
   long  nbvx         = arg(0,stack,0L);
   bool  requireborder= arg(3,stack,arg(1,stack,false));
    KNM<double> * p=0;  p=arg(2,stack,p);
    double alea = arg(4,stack,0.);

   ffassert(   nbvx >= 0);
   return SetAny<pmesh>(Add2StackOfPtr2FreeRC(stack,BuildMesh(stack,borders,false,nbvx,requireborder,p,alea)));

}

AnyType BuildMeshFile::operator()(Stack stack)  const {
     string*  filename = GetAny<string* >((*getfilename)(stack));
    long  nbvx         = arg(0,stack,0);
    ffassert(   nbvx >= 0);
    pmesh pmsh=buildmeshbamg( filename,nbvx);
    Add2StackOfPtr2FreeRC(stack,pmsh);//  07/2008 FH
    return SetAny<pmesh>(pmsh);

}
static int  ChangeLab(const map<int,int> & m,int lab)
{
    map<int,int>::const_iterator i=m.find(lab);
    if(i != m.end())
        lab=i->second;
    return lab;
}

AnyType Op_trunc_mesh::Op::operator()(Stack stack)  const {
    // Remark : F.Hecht feb 2016 ...
    // WARNING for DDM
    // trunc(trunc(Th,op1),op2) =trunc(trunc(Th,op2),op1) => no renumbering ....
    using namespace    Fem2D;
    const pmesh  pTh = GetAny<pmesh>((*getmesh)(stack));
    if( !pTh) return pTh;
    const Mesh & Th = *pTh;
    long kkksplit =std::max(1L, arg(0,stack,1L));
    long label =arg(1,stack,2L);
    KN<long> * pn2o =  arg(2,stack);
    KN<long> * po2n =  arg(3,stack);
    KN<int> split(Th.nt);
    bool renum=arg(4,stack,false);//  change to false too dangerous with ddm the trunc must commute in DDM
    KN<long> kempty;
    KN<long> nre = arg(5,stack,kempty);
    KN<long> nrt = arg(6,stack,kempty);
    Expression flab = nargs[7] ;
    Expression freg = nargs[8] ;
    split=kkksplit;
    long ks=kkksplit*kkksplit;
    MeshPoint *mp= MeshPointStack(stack),mps=*mp;
    long kk=0;
    for (int k=0;k<Th.nt;k++)
     {
       Triangle & K(Th[k]);
       R2 B(1./3.,1./3.);
       mp->set(Th,K(B),B,K,0);
       if (  GetAny<bool>((*bbb)(stack))  ) kk++;
       else  split[k]=0  ;
     }

    if(pn2o)
        {
          pn2o->resize(kk*ks);
          KN<long> &n2o(*pn2o);
          int l=0;
          for(int k=0; k< Th.nt; ++k)
             if( split[k] )
                 for(int i=0; i< ks; ++i)
                     n2o[l++] = k;
        }
        if(po2n)
        {
            po2n->resize(Th.nt);
            KN<long> &o2n(*po2n);
            int l=0;
            for(int k=0; k< Th.nt; ++k)
                if( split[k] )
                {
                        o2n[k] = l;
                       l+=ks;
                }
            else o2n[k]=-1;
        }
//    chnage ....

     *mp=mps;
     if (verbosity>1)
     cout << "  -- Trunc mesh: Nb of Triangle = " << kk << " label=" <<label <<endl;
    Mesh  * pmsh = new Mesh(Th,split,false,label);
    if(renum)
        pmsh->renum();
    {
        Mesh &Th=*pmsh;
        map<int,int> mape,mapt;
        for(int i=0;i<nre.N();i+=2)
            mape[nre[i]]=nre[i+1];
        for(int i=0;i<nrt.N();i+=2)
            mapt[nrt[i]]=nrt[i+1];
        R2 PtHat(1./3,1./3.);
        for (int i=0;i<Th.nt;i++)
        {
            int lab=ChangeLab(mapt,Th[i].lab);
            if(freg)
            {
                mp->set(Th,Th[i](PtHat),PtHat,Th[i],lab);
                Th[i].lab =GetAny<long>( (* freg)(stack)) ;
                 if(verbosity>0) cout << " freg "<<Th[i].lab  << endl;
            }
        }


        // les arete frontieres qui n'ont pas change
        for (int i=0;i<Th.nbBrdElmts();i++)
        {
            int ke,k =Th.BoundaryElement(i,ke);
            int kke,kk= Th.ElementAdj(k,kke=ke);
            const   Triangle &K(Th[k]);
            int l0,l1=ChangeLab(mape,l0=Th.bedges[i].lab) ;
            mp->set(Th,Th[k](PtHat),PtHat,Th[k],l1);
            if(flab)
            {
                R2 E=K.Edge(ke);
                double le = sqrt((E,E));
                double sa=0.5,sb=1-sa;
                R2 PA(TriangleHat[VerticesOfTriangularEdge[ke][0]]),
                PB(TriangleHat[VerticesOfTriangularEdge[ke][1]]);
                R2 Pt(PA*sa+PB*sb ); //
                MeshPointStack(stack)->set(Th,K(Pt),Pt,K,l1,R2(E.y,-E.x)/le,ke);
                Th.bedges[i].lab =GetAny<long>( (*flab)(stack)) ;
            }
        }

    }
    Add2StackOfPtr2FreeRC(stack,pmsh);//  07/2008 FH

  return SetAny<pmesh>(pmsh);
 };


AnyType SplitMesh::operator()(Stack stack) const
{
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
  using  Fem2D::MeshPointStack;
   MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  const  Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
   ffassert(Thh);
   int label=0;
  const  Mesh & Th(*Thh);
   long nbt=Thh->nt;
   KN<int> split(nbt);
   R2 B(1./3.,1./3.);
   int smax=0;
   int smin=100000;
   for (int it=0;it<nbt;it++)
    {
      Triangle & K(Th[it]);
      mp->set(Th,K(B),B,K,K.lab);
      split[it]=GetAny<long>((*U)(stack));
      smin=min(smin,split[it]);
      smax=max(smax,split[it]);
    }
   if(verbosity) cout << "  -- Splitmesh " << Thh << " split  min: " << smin << " max: " << smax << endl;
   Mesh * pth= new Mesh(*Thh,split,false,label);
   R2 Pn,Px;
   pth->BoundingBox(Pn,Px);
   if(!pth->quadtree)
   pth->quadtree=new Fem2D::FQuadTree(pth,Pn,Px,pth->nv);

   *mp=mps;
    Add2StackOfPtr2FreeRC(stack,pth);//  07/2008 FH

    return SetAny<pmesh>(pth);
}
AnyType SaveMesh::operator()(Stack stack) const
{
  using  Fem2D::MeshPointStack;
  const  Fem2D::Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
   string * fn =  GetAny<string*>((*filename)(stack));
   if (!xx && !yy ) {
     const ::bamg::Triangles * bTh= msh2bamg(*Thh);
     (*bTh).Write(fn->c_str(),::bamg::Triangles::AutoMesh);
     delete bTh;
     }
   else {
    MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
     ofstream fp((*fn+".points").c_str());
     ofstream ff((*fn+".faces").c_str());
     fp.precision(12);
     if (verbosity>1)
       cout << "  -- Opening files " << (*fn+".points") << " and " << (*fn+".faces") << endl;
    const   Fem2D::Mesh & Th=*Thh;
    long nbv=Thh->nv;
    long nbt=Thh->nt;
    ff << nbt << endl;
    fp << nbv << endl;
    KN<long> num(nbv);
    num=-1;
    int k=0;
    for (int it=0;it<nbt;it++)
      {
        const Fem2D::Triangle & K(Th[it]);

        num[Th(K[0])]=k++;
        num[Th(K[1])]=k++;
        num[Th(K[2])]=k++;

         ff << " 3 " << Th(K[0])+1 << ' ' << Th(K[1])+1 << ' ' << Th(K[2])+1 << ' '
         	<< " 0 0 0 " << K.lab << '\n';
      }
    if( verbosity>5)
     cout << "  - end writing faces  " << endl;
    for (int iv=0;iv<nbv;iv++)
      {
        const Fem2D::Vertex  & v(Th(iv));
        ffassert( iv == Th(num[iv]/3,num[iv]%3));
        mp->setP(Thh,num[iv]/3,num[iv]%3);

        fp << GetAny<double>((*xx)(stack)) << ' ';
        fp << GetAny<double>((*yy)(stack)) << ' ';
        if (zz)
         fp << GetAny<double>((*zz)(stack)) << ' ';
        else
         fp << " 0 ";
        fp << v.lab<< '\n';
      }
    if( verbosity>5)
     cout << "  - end writing points  " << endl;

   *mp= mps;

   }
  //  delete fn;   modif mars 2006 auto del ptr
   return SetAny<pmesh>(Thh);
}
AnyType MoveMesh::operator()(Stack stack) const
{

  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
  using  Fem2D::MeshPointStack;
   MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
   const Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
   ffassert(Thh);
   long nbv=Thh->nv;
   long nbt=Thh->nt;
   KN<double> u(nbv),v(nbv);
   double infini=DBL_MAX;
   u=infini;
   for (int it=0;it<nbt;it++)
    for (int iv=0;iv<3;iv++)
    {
      int i=(*Thh)(it,iv);
      if ( u[i]==infini) { // if nuset the set
        mp->setP(Thh,it,iv);
        u[i]=GetAny<double>((*U)(stack));
        v[i]=GetAny<double>((*V)(stack));
      }
    }

  const  Mesh * pth= MoveTheMesh(*Thh,u,v);
   if (pth)
     for (size_t i=0;i<sol.size();i++)
       { //  ale
	 pair<FEbase<double,v_fes>,int> * s = GetAny<pair<FEbase<double,v_fes>,int>*>( (*sol[i])(stack));
	 ffassert(s->first.Vh);
	 ffassert( &s->first.Vh->Th == Thh); // same old mesh
	 ffassert(0); // a faire ????
       }
   *mp=mps;
    Add2StackOfPtr2FreeRC(stack,pth);// 07/2008 FH
    return SetAny<pmesh>(pth);

}


AnyType Adaptation::operator()(Stack stack) const
{
  using namespace bamg;
  using bamg::Min;
  using bamg::Max;
  using bamg::Abs;


  Real8 err         = arg(2,stack,0.01);    // coef in the metric
  Real8  errg       = arg(3,stack,Min(0.01,err));// Modif FH 201217
  long nbsx         = Max(100L,arg(4,stack,9000L));
  long nbsmooth     =  arg(5,stack,3L);
  long nbjacobi     =  arg(6,stack,0L) ;              // if increased will be more smooth
  const Real8 raison = arg(7,stack,1.8);
  const Real8 omega =  arg(8,stack,1.0) ;
  bool iso          =   arg(9,stack,false);
  bool AbsError     =   arg(10,stack,true);
  Real8 CutOff      =   arg(11,stack, 1.0e-6);
  verbosity         =   arg(12,stack, (long) verbosity);
  bool inq          =   arg(13,stack,false);
  bool SplitEdgeWith2Boundary =  arg(14,stack,true);
  double maxsubdiv             = Max(Min( arg(15,stack,10.0),10.0),0.1);
  double anisomax              =   Max((double) arg(16,stack,1.0e6),1.0);
  bool rescaling               = arg(17,stack,true) ;
  bool KeepBackVertices        =   arg(18,stack,true) ;
  int givenmetric              =  arg(19,stack,false) ;
  double powerM                = arg(20,stack,1.0) ;
  double cutoffradian          = arg(21,stack,-1.0)* bamg::Pi/180. ;
  bool split                    = arg(22,stack,false) ;
  bool nomeshgeneration         = arg(23,stack,false) ;
  long lwarning=0;

  long * pwarning                = arg(27,stack,&lwarning) ;
  long &warning=*pwarning; //  get get warning message ...

  //   the 24th param is metrix  and is store at compilation time
  //  const E_Array * expmetrix = dynamic_cast<const E_Array *>(nargs[24]);
  //   the 25th param is periodic and it store at compilation time
  // in nbcperiodic,periodic  variable
    KN<long> reqedges0;
    //  list of label of required edges , for no adapattion on this part of the boundary.
    KN<long> reqedges ( nargs[26] ? GetAny< KN_<long> >( (*nargs[26])(stack) ): (KN_<long>)reqedges0);

  if(reqedges.N() && verbosity)
      cout << " reqedges labels "  << reqedges << endl;
  KN<double> *mm11=0, *mm12=0,* mm22=0;

  using Fem2D::MeshPoint;
  using Fem2D::Mesh;
  const  Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
  ffassert(Thh);
   Triangles * oTh =0;
  if (nbcperiodic) {
    KN<int> ndfv(Thh->nv);
    KN<int> ndfe(Thh->neb);
    int nbdfv=0,nbdfe=0;
    BuildPeriodic(nbcperiodic,periodic,*Thh,stack,nbdfv,ndfv,nbdfe,ndfe);
     oTh = msh2bamg(*Thh,cutoffradian,nbdfv,ndfv,nbdfe,ndfe,reqedges,reqedges.N());

  }
  else
   oTh = msh2bamg(*Thh,cutoffradian,reqedges,reqedges.N());
   Triangles &Th(*oTh);
  bool mtx=em11 && em22 && em12;
  if( mtx )
   {
      mm11= GetAny<KN<double> *>( (*em11)(stack) );
      mm22= GetAny<KN<double> *>( (*em22)(stack) );
      mm12= GetAny<KN<double> *>( (*em12)(stack) );
      if (mm11->N() != Th.nbv || mm22->N() != Th.nbv || mm12->N() != Th.nbv)
        ExecError("The size of3  metrics array must be equal to nb of vertex");

   }

  KN<double> &m11=*mm11;
  KN<double> &m12=*mm12;
  KN<double> &m22=*mm22;


  Real8 hmax = 0.3*Th.MaximalHmax(); // final largest edge
  Real8 hmin = Th.MinimalHmin();        // final smallest edge

  Real8 coef =1;                // a priori don't touch
  // gestion des arguments
  hmin              = Max(hmin, arg(0,stack,hmin));
  hmax              = Min(hmax,arg(1,stack,hmax));
  if (inq && initgraphique)
   if (!withrgraphique ) {initgraphique();withrgraphique=true;}

    if (iso)  anisomax=1;
  if (verbosity>2)
    {
      cout << endl  << endl;
      cout << " \t\t ## adapt : nbsol= " << nbsol << ", nbsx = " << nbsx << ", err = " << err ;
      cout << ", hmin = " << hmin << ", hmax = " << hmax <<endl;
      cout << " \t\t    ratio  = " << raison << ", nbsmooth = " << nbsmooth ;
      cout << ", omega = " << omega <<  ", coef = " << coef << ", iso = " << iso << endl;
      cout << " \t\t    AbsError =" << AbsError << ", CutOff = " << CutOff << ", nbjacobi = " << nbjacobi <<endl;
      cout << " \t\t    maxsubdiv = " << maxsubdiv << " splitpbedge = " << SplitEdgeWith2Boundary  <<endl;
      cout << " \t\t    anisomax = " << anisomax << ", rescaling = " << rescaling << ", power = " << powerM
           << ", KeepBackvertices = " << KeepBackVertices << " IsMetric = " << givenmetric
      << endl << endl ;
    }

 Th.ReMakeTriangleContainingTheVertex();
  MeshPoint* mp(Fem2D::MeshPointStack(stack));

  Int4 i,iv;
  int ksol =0;
  for (i=0;i<nbsol;i++)
    ksol += typesol[i]+1; //  marche en 2d

  double * lessol = new double [Th.nbv*ksol];
  double *ss = lessol;
  // be careful because renum --
  // the triangle was no renum
  for ( iv=0;iv<Th.nbv;iv++)
    Th[iv].color=1; // color
 for (Int4  it = 0; it < Thh->nt; it++)
    for (Int4  jt = 0; jt < 3; jt++)
      {
        bamg::Vertex & v= Th(it)[jt];
        if ( v.color)
          {
            v.color =0; // uncolor
            mp->setP(Thh ,it,jt);

            ss = lessol + ksol* Th.Number(v);
            for (int j =0; j < ksol; j++)
              *ss++= GetAny<double>( (*sol[j])(stack) );

          }
      }
  mp->unset();
  // computation of the metric ---
  // better thing -> create keyword in the language
  //    a faire F Hecht .

  Metric Mhmax(hmax);

  for ( iv=0;iv<Th.nbv;iv++)
    Th[iv].m = Mhmax;
    int miss=0;
   if (mtx)
    for ( iv=0;iv<Th.nbv;iv++)
       Th[iv].m.IntersectWith(MetricAnIso(m11[iv],m12[iv],m22[iv]));// add inters ..
      else miss++;
  if ( givenmetric)
    if (ksol == 1)
      {
        for (Int4  iv = 0,k=0; iv < Th.nbv ; iv++)
          Th[iv].m.IntersectWith(Metric(lessol[k++]));
      }
    else if (ksol == 3)
      {
        for (Int4  iv = 0,k=0; iv < Th.nbv ; iv++, k += 3)
          {
	    Metric MM(lessol[k],lessol[k+1],lessol[k+2]);
	    MatVVP2x2 vp(MM);
	    vp.Abs();
	    Th[iv].m.IntersectWith(vp);
          }
      }
    else
      lgerror("Adapt mesh: ksol is wrong, IsMetric  and ksol != 1 or 3");
  else
  Th.IntersectConsMetric(lessol,nbsol,typesol,hmin,hmax,sqrt(err)*coef,anisomax,AbsError?0.0:CutOff,nbjacobi,rescaling,powerM,0);
#ifdef DRAWING1
  if ( (inq!=0) && initgraphique ) {
    if (!withrgraphique) {initgraphique();withrgraphique=true;}
    reffecran();
    Th.InitDraw();
    Th.Draw();
    if (!NoWait) Th.inquire();}
#endif
  delete [] lessol;
  Th.IntersectGeomMetric(errg,iso);
  Th.SmoothMetric(raison);
  Th.MaxSubDivision(maxsubdiv);
  Th.BoundAnisotropy(anisomax);
  // end of metric's computation
   if (mtx)
    for ( iv=0;iv<Th.nbv;iv++)
      {
       m11[iv] = Th[iv].m.a11  ;
       m22[iv] = Th[iv].m.a22 ;
       m12[iv] = Th[iv].m.a21;
      }

  Triangles* nTh = 0;

  if ( ! nomeshgeneration)
  {
  nTh= new Triangles(nbsx,Th,KeepBackVertices); // Adaption is here
  if (split)
	    nTh->SplitElement(1); // modif FH mai 2009 (thank J-M Mirebeau) : Th ->nTh

  if(SplitEdgeWith2Boundary)
    nTh->SplitInternalEdgeWithBorderVertices();
  if(verbosity>3)
    nTh->ShowHistogram();
  if (nbsmooth)
    nTh->SmoothingVertex(nbsmooth,omega);
  if(verbosity>2 && nbsmooth)
    nTh->ShowHistogram();
  if(verbosity>0)
      nTh->ShowRegulaty()  ;

#ifdef DRAWING
  if ((inq!=0) && initgraphique) {
    if (!withrgraphique ) {initgraphique();withrgraphique=true;}

    reffecran();
    nTh->InitDraw();
    nTh->Draw();
    if(!NoWait) nTh->inquire();}
#else
  inq=0;
#endif
  Metric M(hmax);
  for (iv=0;iv < Th.nbv;iv++)
    Th[iv].m = M;
  warning = nTh->warning;

  const Mesh * g=  bamg2msh(nTh,true);

  delete nTh;
  delete oTh;
  Add2StackOfPtr2FreeRC(stack,g);// 07/2008  FH

  return SetAny<pmesh>(g);
  }
 else {

     if(verbosity>1)
     {  cout << " regularty Old mesh / New metrix ";
	 oTh->ShowRegulaty()  ;}
   delete oTh;
   return SetAny<pmesh>(Thh);
 }
}

const Fem2D::Mesh  * EmptyTheMesh(const Fem2D::Mesh *  const & pTh,long *ssd=0)
{
  using namespace Fem2D;
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
  const Mesh & Th=*pTh;
  using  Fem2D::MeshPointStack;
  int nbv=Th.nv;
  int nbt=0;
  int neb=Th.neb;
  KN<int> renum(Th.nv);
  renum=-1;
  if(ssd)
   {
    nbt=-1; //  pour ne pas retire l'exterieur
    int nebb=0;
    for (int i=0;i<Th.nt;i++)
     for (int e=0;e<3;e++)
      { int ee=e,ii=Th.ElementAdj(i,ee);
        if ( (ii>= 0 && ii != i) && (i<ii && ssd[ii] != ssd[i]) )
          {
           nebb++;
           int i1= Th(i,VerticesOfTriangularEdge[e][0]);
           int i2= Th(i,VerticesOfTriangularEdge[e][1]);
           cout << i1 << " " << i2 << endl;
           renum[i1]=1;
           renum[i2]=1;
         }
      }
      neb=nebb;
   }
   else
  for (int i=0;i<neb;i++)
    {
     int i1=Th(Th.bedges[i][0]);
     int i2=Th(Th.bedges[i][1]);
     renum[i1]=1;
     renum[i2]=1;
     }
  int nbvnew=0;
  for(int i=0;i<nbv;i++)
   if (renum[i]>=0)
     renum[i]= nbvnew++;

  Vertex * v= new Vertex[nbvnew];
  //  Triangle *t= 0;
  BoundaryEdge *b= new BoundaryEdge[neb];
  //Vertex *vv=v;
  Vertex *vo=Th.vertices;
  BoundaryEdge * bb=b;

  if(ssd)
   {
    int nebb=0;
    for (int i=0;i<Th.nt;i++)
     for (int e=0;e<3;e++)
      { int ee=e,ii=Th.ElementAdj(i,ee);
        if ( (ii>= 0 && ii != i) && (i<ii && ssd[ii] != ssd[i]) )
         {
           nebb++;
           int i1= renum[Th(i,VerticesOfTriangularEdge[e][0])];
           int i2= renum[Th(i,VerticesOfTriangularEdge[e][1])];
           int labM=0,labm=0;
           labM=Th[i].lab;
           if( ii >=0) labm=Th[ii].lab;
           if( labM <labm) Exchange(labM,labm);
           int lab= 100*labM+labm ;
           *bb++ = BoundaryEdge(v,i1,i2,lab);
         }
      }
   }
   else
  for (int i=0;i<neb;i++)
    {
     int i1=renum[Th(Th.bedges[i][0])];
     int i2=renum[Th(Th.bedges[i][1])];
     int lab=Th.bedges[i].lab;
     *bb++ = BoundaryEdge(v,i1,i2,lab);
    }

  for (int i=0;i<nbv;i++)
   {
     int j=renum[i];
     if(j>=0)
      {
        v[j].x=vo[i].x;
        v[j].y=vo[i].y;
        v[j].lab=vo[i].lab;
      }
    }

 {
  Mesh * m = new Mesh(nbvnew,nbt,neb,v,0,b);
  R2 Pn,Px;
  m->BoundingBox(Pn,Px);
  m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv);
//  m->decrement();
  return m;

}
}
const Fem2D::Mesh  * EmptyTheMesh( const Fem2D::Mesh *  const & pTh)
{
  return EmptyTheMesh(pTh,0);
}
const Fem2D::Mesh  * EmptyTheMesh(const  Fem2D::Mesh *  const & pTh,  KN<long> * const &  k)
{
  return EmptyTheMesh(pTh,*k);
}


const Mesh * MoveTheMesh(const Fem2D::Mesh &Th,const KN_<double> & U,const KN_<double> &V)
{
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
 // using  Fem2D::R;
  using  Fem2D::MeshPointStack;
  int nbv=Th.nv;
  int nbt=Th.nt;
  int neb=Th.neb;
  Vertex * v= new Vertex[nbv];
  Triangle *t= new Triangle[nbt];
  BoundaryEdge *b= new BoundaryEdge[neb];
  Vertex *vv=v;
  for (int i=0;i<nbv;i++)
   {
     R2 P(U[i],V[i]);
     vv->x=P.x;
     vv->y=P.y;
     vv->lab = Th(i).lab;
     vv++;
    }
  Triangle *tt= t;
  int nberr=0;
  double atotal=0,atotal0=0;
  double amax=-1e100,amin=1e100;
  for (int i=0;i<nbt;i++)
    {
       atotal0 += Th[i].area;
      int i0=Th(i,0), i1=Th(i,1),i2=Th(i,2);
      double a= Area2(v[i0],v[i1],v[i2])/2;
      atotal +=  a;
      amax=max(a,amax);
      amin=min(a,amin);
    }
  double  eps = 1e-6 * max(abs(atotal),1e-100)/atotal0;
  bool rev=(atotal<0);
  if(verbosity>2 && rev) cout <<  "  -- movemesh negatif tranfomation => reverse all triangle (old area " << atotal0 << ") ( new area  " << atotal << ") "<< endl;
  for (int i=0;i<nbt;i++)
    {
      int i0=Th(i,0), i1=Th(i,1),i2=Th(i,2);
      if(rev) swap(i1,i2);
      R a= Area2(v[i0],v[i1],v[i2])/2;
      if ( a < Th[i].area*eps)
       { nberr++;
        if (verbosity>1)
         {
          if (nberr==1) { cerr << "Erreur: MoveMesh "; }
	  cerr << " T = " << Th[i] <<  endl;
          }
          if (nberr < verbosity*5) {
            cerr << " " <<i;
            if ( nberr % 5 )  cerr << "\n\t";}
         }
       else
        (*tt++).set(v,i0,i1,i2,Th[i].lab,a);

   }
    if (nberr)
      { if (verbosity)
	  cerr << "Error movemesh: " << nberr << " triangles was reverse  (=> no move)" <<  endl;
	  cout << " u min " << U.min() << " max " << U.max() << endl;
	  cout << " v min " << V.min() << " max " << V.max() << endl;

	  delete []v;
	  delete []t;
	  delete []b;
	  throw(ErrorExec("Error move mesh triangles was reverse",1));
	  return 0;
      }

  BoundaryEdge * bb=b;
  for (int i=0;i<neb;i++)
    {
     int i1=Th(Th.bedges[i][0]);
     int i2=Th(Th.bedges[i][1]);
     int lab=Th.bedges[i].lab;
     if(rev) swap(i1,i2);
     *bb++ = BoundaryEdge(v,i1,i2,lab);
    }
 {
  Mesh * m = new Mesh(nbv,nbt,neb,v,t,b);
  R2 Pn,Px;
  m->BoundingBox(Pn,Px);
  m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv);
  return m;

}
}

/// <<Carre>> Builds a square-shaped 2D mesh. An Expression [[file:AFunction.hpp::Expression]] is a pointer to an object
/// of class E_F0 [[file:AFunction.hpp::E_F0]].
Mesh * Carre(int nx,int ny,Expression fx,Expression fy,Stack stack,int flags,KN_<long> lab,long reg)
{
    Mesh * m=Carre_( nx, ny, fx, fy, stack, flags,lab, reg);
    Add2StackOfPtr2FreeRC(stack,m);// 07/2008 FH
    return m;

}
Mesh * Carre_(int nx,int ny,Expression fx,Expression fy,Stack stack,int flags,KN_<long> lab,long reg)
{
  if(verbosity>99)  cout << " region = " << reg << " labels " << lab <<endl;
  const int unionjack=1;
  using  Fem2D::Vertex;
  using  Fem2D::Triangle;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
  using  Fem2D::R;

  using  Fem2D::MeshPointStack;
  int nx1=nx+1,ny1=ny+1;
  int nbv=(nx1)*(ny1);
  int nbt=(nx)*ny*2;
  int neb=(nx+ny)*2;
  int l1=1,l2=2,l3=3,l4=4;
  if(lab.N()>=1) l1=lab[0];
  if(lab.N()>=2) l2=lab[1];
  if(lab.N()>=3) l3=lab[2];
  if(lab.N()>=4) l4=lab[3];

  Vertex * v= new Vertex[nbv];
  Triangle *t= new Triangle[nbt];
  BoundaryEdge *b= new BoundaryEdge[neb];
  Vertex *vv=v;
  for (int j=0;j<ny1;j++)
   for (int i=0;i<nx1;i++)
   {
     R2 P((R) i/nx,(R) j/ny);

     vv->x=P.x;
     vv->y=P.y;
     vv->lab =  0;
     vv++;
   }
  if (fx || fy)
   {
     MeshPoint *mp=MeshPointStack(stack);
    vv=v;
    for (int j=0;j<ny1;j++)
     for (int i=0;i<nx1;i++)
       {
         mp->set(vv->x,vv->y);
         if (fx)  vv->x= GetAny<R>((*fx)(stack));
         if (fy)  vv->y= GetAny<R>((*fy)(stack));
         vv++;
       }
   }
  Triangle *tt= t;
  //
  bool direct = det(v[0],v[1],v[nx1+1]) > 0; //  signe  triangle 0

    if(verbosity>1&& !direct) cout << "  -- square : all triangles are reversed" <<   endl;

    int p[2]={1,0};
    if(direct) {p[0]=0,p[1]=1;}
  for (int j=0;j<ny;j++)
   for (int i=0;i<nx;i++)
     {
       int i0 = i + j*nx1;
       int i1= i0+1;
       int i2=i1+nx1;
       int i3=i2-1;
       bool c00= (i==0) && (j==0);
       bool c10= (i==nx-1) && (j==0);
       bool c01= (i==0) && (j==ny-1);
       bool c11= (i==nx-1) && (j==ny-1);

      bool cas = true;
      switch (flags)
       {
       case unionjack:
          cas =  (i+ j) %2 ; break;
       case 2: // new  jan 2010
	     cas = false; break;
       case 3:// never triangle with 3 vertex on boundary
	     cas = c01 || c10 ? false : true; break;
       case 4:// never triangle with 3 vertex on boundary
	     cas = c00 || c11 ? true : false; break;
       default:
          cas = true;
       } ;

	 if (cas)
	   {
	       if(direct)
		 {     // diag 1
		     (tt++)->set(v,i0,i1,i2,reg,0.0);
		     (tt++)->set(v,i0,i2,i3,reg,0.0);
		 }
	       else
		 {     // diag 1
		     (tt++)->set(v,i1,i0,i2,reg,0.0);
		     (tt++)->set(v,i2,i0,i3,reg,0.0);
		 }
	   }
	 else
	   {
	       if(direct)
		 {    // diag 2
		     (tt++)->set(v,i0,i1,i3,reg,0.0);
		     (tt++)->set(v,i3,i1,i2,reg,0.0);
		 }
	       else
		 {    // diag 2
		     (tt++)->set(v,i1,i0,i3,reg,0.0);
		     (tt++)->set(v,i1,i3,i2,reg,0.0);

		 }
	   }

     }
    BoundaryEdge * bb=b;
  for (int i=0;i<nx;i++)
    {  // bottom
      int i1=i,i2=i1+1;
      *bb++ = BoundaryEdge(v,i1,i2,l1);
      v[i1].lab=v[i2].lab=l1;
    }
  for (int j=0;j<ny;j++)
    { // right
      int i=nx;
      int i1= i + j*nx1 ,i2=i1+nx1;
      *bb++ = BoundaryEdge(v,i1,i2,l2);
      v[i1].lab=v[i2].lab=l2;
    }

  for (int i=0;i<nx;i++)
    {  // up
      int j=ny;
      int i1=i + j*nx1,i2=i1+1;
      *bb++ = BoundaryEdge(v,i1,i2,l3);
      v[i1].lab=v[i2].lab=l3;

    }
  for (int j=0;j<ny;j++)
    { // left
      int i=0;
      int i1= i + j*nx1,i2=i1+nx1;
      *bb++ = BoundaryEdge(v,i1,i2,l4);
      v[i1].lab=v[i2].lab=l4;
    }

  {
    if(verbosity) cout << "  -- Square mesh : nb vertices  =" << nbv
		       << " ,  nb triangles = " << nbt << " ,  nb boundary edges " << neb << endl;
    Mesh * m = new Mesh(nbv,nbt,neb,v,t,b);
    R2 Pn,Px;
    m->BoundingBox(Pn,Px);
    m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv);
    return m;
  }
}

/// <<MeshCarre2>> Creates a square mesh by calling the [[Carre]] function at script evaluation time. Uses class E_F0mps
/// [[file:AFunction.hpp::E_F0mps]] to connect to the FF language.

class MeshCarre2 :   public E_F0mps {
public:
  typedef pmesh  Result;
  Expression nx,ny;
  static basicAC_F0::name_and_type name_param[] ;
  static const int n_name_param =1+2;
  Expression nargs[n_name_param];

  long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}


  MeshCarre2(const basicAC_F0 & args)
  {
    args.SetNameParam(n_name_param,name_param,nargs);
    nx=to<long>(args[0]);
    ny=to<long>(args[1]);
  }

  static ArrayOfaType  typeargs() {
    return  ArrayOfaType(atype<long>(),atype<long>(),false);}

  /// <<MeshCarre2_f>>

  static  E_F0 * f(const basicAC_F0 & args){
    return new MeshCarre2(args);}

  AnyType operator()(Stack s) const {
    long flags=arg(0,s,0);
    KN<long> zz;
    KN<long> label=arg(1,s,zz);
    long region=arg(2,s,0L);// correct aout 2010 FH ... 2-> 0

    /// calls [[Carre]]

    return SetAny<pmesh>(Carre( GetAny<long>( (*nx)(s)) , GetAny<long>( (*ny)(s)),0,0,s,flags,label,region ));
  }

  operator aType () const { return atype<pmesh>();}
};

/// <<MeshCarre2f>> Creates a square mesh by calling the [[Carre]] function at script evaluation time

class MeshCarre2f :   public E_F0mps {
public:
  typedef pmesh  Result;
  Expression nx,ny;
  Expression fx,fy;
  static basicAC_F0::name_and_type name_param[] ;
  static const int n_name_param =1+2;
  Expression nargs[n_name_param];

  long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}

  MeshCarre2f(const basicAC_F0 & args)
  {
    args.SetNameParam(n_name_param,name_param,nargs);
    nx=to<long>(args[0]);
    ny=to<long>(args[1]);
    const E_Array *  a= dynamic_cast<const E_Array*>(args[2].LeftValue());
    ffassert(a);fx=0;fy=0;
    if (a->size()>0) fx=to<double>( (*a)[0]);
    if (a->size()>1) fy=to<double>( (*a)[1]);
  }

  static ArrayOfaType  typeargs() {
    return  ArrayOfaType(atype<long>(),atype<long>(),atype<E_Array>(),false);}

  static  E_F0 * f(const basicAC_F0 & args){
    return new MeshCarre2f(args);}

  AnyType operator()(Stack s) const {
    long flags=arg(0,s,0);
    KN<long> zz;
    KN<long> label=arg(1,s,zz);
    long region=arg(2,s,0L);

    return SetAny<pmesh>(Carre( GetAny<long>( (*nx)(s)) , GetAny<long>( (*ny)(s)), fx,fy,s , flags,label,region));}

  operator aType () const { return atype<pmesh>();}
};

basicAC_F0::name_and_type  MeshCarre2::name_param[]= {
	{  "flags", &typeid(long) },
	{  "label", &typeid(KN_<long> ) },
	{  "region", &typeid(long) }

};
basicAC_F0::name_and_type  MeshCarre2f::name_param[]= {
	{  "flags", &typeid(long) },
	{  "label", &typeid(KN_<long> )},
	{  "region", &typeid(long)}
    };

   void  MeshErrorIO(ios& )
{
   ExecError("Mesh IO Error ");
}

inline pmesh *  initMesh(pmesh * const & p, string * const & s) {
  Mesh * m;
  *p= m =new Mesh(*s);
  m->MakeQuadTree();
  return p;
 }

class CheckMoveMesh :  public E_F0mps { public:

   typedef double  Result;
   Expression getmesh;
   Expression U,V;
   int nbsol;
    vector<Expression> sol;

    CheckMoveMesh(const basicAC_F0 & args) :nbsol(args.size()-2),sol(args.size()-2)
    {
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
  using  Fem2D::MeshPointStack;
      args.SetNameParam();
      getmesh=to<pmesh>(args[0]);
      const E_Array * a = dynamic_cast<const E_Array *>(args[1].LeftValue());

      ffassert(a);
      if (a->size() !=2) CompileError("CheckMoveMesh(Th,[u,v]) need 2 componate in array ",atype<pmesh>());
      U=to<double>( (*a)[0]);
      V=to<double>( (*a)[1]);

      for (int i=2;i<args.size();i++)
        sol[i-2]=to<double>(args[i]);
    }
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),atype<E_Array>(),false);}
    static  E_F0 * f(const basicAC_F0 & args){ return new CheckMoveMesh(args);}
    AnyType operator()(Stack s) const ;
    operator aType () const { return atype<double>();}

};
AnyType CheckMoveMesh::operator()(Stack stack) const
{

  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
  using  Fem2D::MeshPointStack;
   MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
   const Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
   const Mesh & Th(*Thh);
   ffassert(Thh);
   long nbv=Thh->nv;
   long nbt=Thh->nt;
   KN<double> u(nbv),v(nbv);
   double infini=DBL_MAX;
   u=infini;
   for (int it=0;it<nbt;it++)
    for (int iv=0;iv<3;iv++)
    {
      int i=(*Thh)(it,iv);
      if ( u[i]==infini) { // if nuset the set
        mp->setP(Thh,it,iv);
        u[i]=GetAny<double>((*U)(stack));
        v[i]=GetAny<double>((*V)(stack));
      }
    }
     double minarea=DBL_MAX;
    for (int t=0;t<Th.nt;t++)
     {
      int i0=Th(t,0),i1=Th(t,1),i2=Th(t,2);
      minarea=Min(minarea,Area2(R2(u[i0],v[i0]), R2(u[i1],v[i1]),R2(u[i2],v[i2])));
     }
    *mp=mps;
    return SetAny<double>(minarea/2.);

}

bool SameMesh(const Mesh * const & pTh1,const Mesh * const & pTh2)
{
    typedef Mesh::Element Element;
    if( !pTh1) return 0;
    if( !pTh2) return 0;
    if( pTh1 == pTh2) return 1;
    if( pTh1->nv != pTh2->nv) return 0;
    if( pTh1->nt != pTh2->nt) return 0;
  //  const Mesh & Th1=*pTh1, & Th2 = *pTh2;
    ffassert(0); // a faire..

    return 1;
}

bool AddLayers(Mesh const * const & pTh, KN<double> * const & psupp, long const & nlayer,KN<double> * const & pphi)
{
    ffassert(pTh && psupp && pphi);
    const int nve = Mesh::Element::NbV;
    const Mesh & Th= *pTh;
    const int nt = Th.nt;
    const int nv = Th.nv;

    KN<double> & supp(*psupp);
    KN<double> u(nv), s(nt);
    KN<double> & phi(*pphi);
    ffassert(supp.N()==nt);//P0
    ffassert(phi.N()==nv); // P1
    s = supp;
    phi=0.;

    for(int step=0; step < nlayer; ++ step)
    {
        u = 0.;
        for(int k=0; k<nt; ++k)
            if(s[k] > 0.0)
                for(int i=0; i<nve; ++i)
                    u[Th(k,i)] = 1.0;

        phi += u;

        s = 0.;
        for(int k=0; k<nt; ++k)
            for(int i=0; i<nve; ++i)
                if(u[Th(k,i)] > 0.0)
                    s[k] = 1.0;

        supp += s;
    }
    phi *= (1./nlayer);
    return true;
}


double arealevelset(Mesh const * const & pTh,KN<double>  * const & pphi,const double & phi0,KN<double>  * const & where)
{
    const Mesh & Th=*pTh;
    KN<double> & phi=*pphi;
    ffassert( phi.N() == Th.nv); // P1
    double arean=0., areap=0.;
    for (int k=0;k<Th.nt;k++)
    {
        if( !where || (*where)[k] )
        {
            const Triangle & K(Th[k]);
            int iK[3]={Th(k,0),Th(k,1),Th(k,2)};
            R  fk[3]={phi[iK[0]]-phi0,phi[iK[1]]-phi0,phi[iK[2]]-phi0 };
            int i0 = 0, i1 = 1, i2 =2;
            if( fk[i0] > fk[i1] ) swap(i0,i1) ;
            if( fk[i0] > fk[i2] ) swap(i0,i2) ;
            if( fk[i1] > fk[i2] ) swap(i1,i2) ;

            if ( fk[i2] <=0.) arean+= K.area;
            else  if ( fk[i0] >=0.) areap+= K.area;
            else {
                double c = (fk[i2]-fk[i1])/(fk[i2]-fk[i0]); // coef Up Traing
                if( fk[i1] < 0 ) { double y=fk[i2]/(fk[i2]-fk[i1]); c *=y*y; }
                else {double y=fk[i0]/(fk[i0]-fk[i1]) ; c = 1.- (1.-c)*y*y; };
                assert( c > -1e-10 && c-1. < 1e-10);
                areap += c*K.area;
                arean += (1-c)*K.area;
            }
        }
    }
    return arean; //  negative area ...
}
double arealevelset(Mesh const * const & pTh,KN<double>  * const & pphi,const double & phi0)
{
    return arealevelset(pTh,pphi,phi0,0);
}


double volumelevelset(Mesh3 const * const & pTh,KN<double> * const &pphi,const double & phi0,KN<double> * const & where)
{
    const Mesh3 & Th = *pTh;
    KN<double> & phi = *pphi;
    ffassert(phi.N() == Th.nv); //assertion fails if the levelset function is not P1
    double volumen=0.,volumep=0.;
    for(int k=0;k<Th.nt;++k)
    {
        if(!where || (*where)[k])
        {
            const Tet & K(Th[k]);
            int iK[4] =  {Th(k,0),Th(k,1),Th(k,2),Th(k,3)};
            R  fk[4]={phi[iK[0]]-phi0,phi[iK[1]]-phi0,phi[iK[2]]-phi0,phi[iK[3]]-phi0};
            int ij[4] = {0,1,2,3};
            for(int i=0;i<4;++i) for(int j=i+1;j<4;++j) if(fk[ij[i]] > fk[ij[j]]) swap(ij[i],ij[j]);
            if(fk[ij[3]]<=0.) volumen += K.mesure();
            else if(fk[ij[0]]>=0.) volumep += K.mesure();
            else
            {
                const R3 A[4] = {K[ij[0]],K[ij[1]],K[ij[2]],K[ij[3]]};
                const R f[4] = {fk[ij[0]],fk[ij[1]],fk[ij[2]],fk[ij[3]]};
                if(f[0]<0. && f[1]>=0.) //the only negative dof is f0
                {
                    double c = (f[0]*f[0]*f[0])/((f[0]-f[1])*(f[0]-f[2])*(f[0]-f[3]));
                    volumen += c*K.mesure();
                    volumep += (1.-c)*K.mesure();
                }
                else if(f[1]<0. && f[2]>=0.) //two negative dof, this is the hard case...
                {
                    const R3  A03 = A[3] + f[3]/(f[3]-f[0]) * (A[0]-A[3]) , A13 = A[3] + f[3]/(f[3]-f[1]) * (A[1]-A[3]),
                    A02 = A[2] + f[2]/(f[2]-f[0]) * (A[0]-A[2]) , A12 = A[2] + f[2]/(f[2]-f[1]) * (A[1]-A[2]);
                    double V1 = f[3]*f[3]/((f[3]-f[1])*(f[3]-f[0])) * K.mesure();
                    double V2 = 1./6.*abs(det(A13-A[2],A12-A[2],A02-A[2]));
                    double V3 = 1./6.*abs(det(A13-A[2],A02-A[2],A03-A[2]));
                    volumep += V1+V2+V3;
                    volumen += (K.mesure() - V1-V2-V3);
                }
                else
                {
                    double c = (f[3]*f[3]*f[3])/((f[3]-f[0])*(f[3]-f[1])*(f[3]-f[2]));
                    volumen += (1.-c)*K.mesure();
                    volumep += c*K.mesure();
                }
            }
        }
    }
    return volumen;
}
double volumelevelset(Mesh3 const * const &pTh,KN<double> * const &pphi,const double & phi0) {return volumelevelset(pTh,pphi,phi0,0);}

long Boundingbox(KN<double>* const& pb, pmesh const& pTh)
{
     KN<double> & bb =*pb;
     if(pTh && bb.N()>=4)
     {
         R2 Pn,Px;
         pTh->BoundingBox(Pn,Px);
         bb[0] = Pn.x;
         bb[1] = Px.x;
         bb[2] = Pn.y;
         bb[3] = Px.y;
         return 0;
     }
    return -1; // error
}
template < class ppmesh >
long Boundingbox(KN<double>* const& pb, ppmesh const& pTh)
{
    KN<double> & bb =*pb;
    if(pTh && bb.N()>=6)
    {
        R3 Pn=pTh->Pmin,Px=pTh->Pmax  ;
        bb[0] = Pn.x;
        bb[1] = Px.x;
        bb[2] = Pn.y;
        bb[3] = Px.y;
        bb[4] = Pn.z;
        bb[5] = Px.z;
        return 0;
     }
    return -1; // error
}

long Boundingbox(pmeshL const& pTh,KN<double>* const& pb )
{ return  Boundingbox<pmeshL>(pb,pTh);}

long Boundingbox(pmeshS const& pTh,KN<double>* const& pb )
{ return  Boundingbox<pmeshS>(pb,pTh);}

long Boundingbox(pmesh3 const& pTh,KN<double>* const& pb )
{ return  Boundingbox<pmesh3>(pb,pTh);}

long Boundingbox(pmesh const& pTh,KN<double>* const& pb )
{ return  Boundingbox(pb,pTh);}
double Chi(Stack stack,pmesh const &pTh)
{  // version 2d  oct 2017 FH.
    if(pTh == 0) return 0.;
    R2 PHat;
    bool outside;
    MeshPoint & mp = *MeshPointStack(stack);
    if(pTh == mp.Th) return 1.;// point of mesh
    const Triangle * K=pTh->Find(mp.P.p2(),PHat,outside);
    if (!outside)
        mp.set(*pTh,mp.P.p2(),PHat,*K,K->lab);
    else return 0.;
    return 1.;
}
double Chi(Stack stack,pmesh3 const &pTh)
{ // version 3d oct 2017 FH.
    if(pTh == 0) return 0.;
    R3 PHat;
    bool outside;

    MeshPoint & mp = *MeshPointStack(stack);

    if(pTh == mp.Th3) return 1.;
    const Tet * K=pTh->Find(mp.P,PHat,outside);
    if (!outside)
        mp.set(*pTh,mp.P,PHat,*K,K->lab);
    else return 0.;
    return 1.;
}
long savegnuplot(pmesh pTh,string* pgp)
    {
        const  Mesh &Th(*pTh);
        const string &gp=*pgp;
        {
            ofstream of(gp.c_str());
            for(int k=0; k<Th.nt;++k)
            {
                R2 G= ((R2) Th[k][0] + Th[k][1] + Th[k][2])/3.;
                of << G << " "<< k << "\n\n\n";
                for(int ip=0; ip<=3;++ip)
                {
                    int i3=ip%3;
                    int i = Th(k,i3);
                    of << (R2) Th[k][i3] << " " << i << endl;
                }
                of << "\n\n";
            }
        }
        if(verbosity>1)
        {
        cout << " to plot with gnuplot plot do under  gnuplot " <<endl;
        cout << " plot '"<<gp<<"' w l,'' w labels offset 1."<< endl;
        }
        return 0;
    }
extern void init_glumesh2D();

void init_lgmesh() {
  if(verbosity&&(mpirank==0) )  cout <<"lg_mesh ";
  bamg::MeshIstreamErrorHandler = MeshErrorIO;
  Global.Add("buildmesh","(",new OneOperatorCode<classBuildMesh>);
  Global.Add("buildmesh","(",new OneOperatorCode<classBuildMeshArray>);
  Global.Add("buildmesh","(",new OneOperatorCode<BuildMeshFile>);

  Global.Add("buildmeshborder","(",new OneOperator1s_<pmesh,const E_BorderN *>(BuildMeshBorder));
  Global.Add("adaptmesh","(",new OneOperatorCode<Adaptation>);
  Global.Add("movemesh","(",new OneOperatorCode<MoveMesh>);
  Global.Add("splitmesh","(",new OneOperatorCode<SplitMesh>);
  Global.Add("checkmovemesh","(",new OneOperatorCode<CheckMoveMesh>);

  /// <<square_keyword>> see [[file:AFunction.hpp::OneOperatorCode]]
  Global.Add("square","(",new OneOperatorCode<MeshCarre2>);
  Global.Add("square","(",new OneOperatorCode<MeshCarre2f>);

  Global.Add("savemesh","(",new OneOperatorCode<SaveMesh>);
  Global.Add("trunc","(", new Op_trunc_mesh);
  Global.Add("readmesh","(",new OneOperator1_<pmesh,string*, E_F_F0_Add2RC<pmesh,string*> >(ReadMeshbamg));
  Global.Add("emptymesh","(",new OneOperator1_<pmesh,pmesh, E_F_F0_Add2RC<pmesh,pmesh> >(EmptyTheMesh));
  Global.Add("emptymesh","(",new OneOperator2_<pmesh,pmesh,KN<long> *, E_F_F0F0_Add2RC<pmesh,pmesh,KN<long>*> >(EmptyTheMesh));
  Global.Add("triangulate","(",new OneOperator1_<pmesh,string*, E_F_F0_Add2RC<pmesh,string*> >(ReadTriangulate));
  Global.Add("triangulate","(",new OneOperator2_<pmesh,KN_<double>,KN_<double>,E_F_F0F0_Add2RC<pmesh,KN_<double>,KN_<double>,E_F0> >(Triangulate));
  TheOperators->Add("<-",
		    new OneOperator2_<pmesh*,pmesh*,string* >(&initMesh));
    // Thg,suppi[],nnn,unssd[]
  Global.Add("AddLayers","(",new OneOperator4_<bool,const Mesh * , KN<double> * , long ,KN<double> * >(AddLayers));
  Global.Add("SameMesh","(",new OneOperator2_<bool,const Mesh * ,const  Mesh * >(SameMesh));
  // use for :   mesh Th = readmesh ( ...);
    //  Add FH mars 2015 to compute mesure under levelset ...
    Global.Add("arealevelset","(",new OneOperator3_<double,const Mesh *, KN<double>*,double>(arealevelset));
    Global.Add("arealevelset","(",new OneOperator4_<double,const Mesh *, KN<double>*,double,KN<double>*>(arealevelset));
    Global.Add("volumelevelset","(",new OneOperator3_<double,const Mesh3 *,KN<double>*,double>(volumelevelset));
    Global.Add("volumelevelset","(",new OneOperator4_<double,const Mesh3*,KN<double>*,double,KN<double>*>(volumelevelset));
    // add FH to get bounding box ...
    Global.Add("boundingbox", "(", new OneOperator2_<long, KN<double>*, pmesh>(Boundingbox));
    Global.Add("boundingbox", "(", new OneOperator2_<long, KN<double>*, pmesh3>(Boundingbox));
    Global.Add("boundingbox", "(", new OneOperator2_<long, KN<double>*, pmeshS>(Boundingbox));
    Global.Add("boundingbox", "(", new OneOperator2_<long, pmesh,KN<double>*>(Boundingbox));
    Global.Add("boundingbox", "(", new OneOperator2_<long, pmesh3,KN<double>*>(Boundingbox));
    Global.Add("boundingbox", "(", new OneOperator2_<long, pmeshS,KN<double>*>(Boundingbox));
    Global.Add("boundingbox", "(", new OneOperator2_<long, pmeshL,KN<double>*>(Boundingbox));
    Global.Add("chi", "(", new OneOperator1s_<double,pmesh,E_F_F0s_<double,pmesh,E_F0mps> >(Chi)); // oct 2017 FH function characteristic
    Global.Add("chi", "(", new OneOperator1s_<double,pmesh3,E_F_F0s_<double,pmesh3,E_F0mps> >(Chi));// oct 2017 FH function characteristic
    Global.Add("savegnuplot","(",new OneOperator2<long,pmesh,string*>(savegnuplot));

  TheOperators->Add("<-",
		    new OneOperator2_<pmesh*,pmesh*,pmesh >(&set_copy_incr));
  init_glumesh2D();
}
