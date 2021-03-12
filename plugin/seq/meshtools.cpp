/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

// to compile: ff-c++ bstream.cpp
// WARNING: do not compile under windows

#include "ff++.hpp"
#include "AFunction_ext.hpp"

template<class Mesh> struct nbvertexElement { static const  int value = Mesh::Element::nv; };
template<> struct nbvertexElement<Mesh> { static const  int value = Mesh::Element::NbV; };

template<class Mesh,class K>
long connexecomponantea(const Mesh * pTh,KN<K>* pcc)
{
    const Mesh & Th= *pTh;
    const long nv = Th.nv;
    const long nt = Th.nt;
    const int d = Mesh::Rd::d;
    const int dh = Mesh::RdHat::d;

    const long  nvk = nbvertexElement<Mesh>::value;
   if(verbosity>9)  cout << " nvk ="<< nvk << endl;
    if(pcc->N() != Th.nt) pcc->resize(Th.nt);
    KN<K> & cc=*pcc;
    long nbc =Th.nt;
    KN<long> forest(Th.nt,-1);// KrusKal algo
    for( int k=0; k< Th.nt;++k)
      for( int e =0; e< nvk; ++e )
       {
           // element adj
           int ee=e,kk = Th.ElementAdj(k,ee);
           if( kk== k || kk<0) continue; // no ajd .
           //  k and kk are in the same cnx. comp.
           long r1=k,r2=kk;
           // recherche des racines
           while( forest[r1]>=0) r1 = forest[r1];
           while( forest[r2]>=0) r2 = forest[r2];
           if(r1 !=r2)
           {
               nbc--;
               if(forest[r1]<forest[r2])
                   forest[r2]=r1; // lie  arbre r2 sur racine r1
               else if(forest[r2]<forest[r1])
                   forest[r1]=r2; // lie  arbre r2 sur racine r1
               else
                   forest[r1]=r2, forest[r2]--;  // lie  arbre r2 sur racine r1 et incremente al profondeur
           }
        }
    //  build cc
    cc=-1; //
    long nc =0;
    for(long s=0; s< nt; ++s)
    {
        long r=s;
        while( forest[r]>=0) r = forest[r];
        if( cc[r] <0) cc[r] = nc++;
        cc[s]=cc[r];
    }
    ffassert(nc==nbc);
    if( verbosity) cout << "  The number of  connexe componante (by adj)  Mesh "<< pTh << " is " << nc << " / "
        " dim = " << d << " dim s " << dh <<  endl;
    return nc;
}
template<class Mesh,class K>
long connexecomponantev(const Mesh * pTh,KN<K>* pcc)
{
    const Mesh & Th= *pTh;
    const long nv = Th.nv;
    const long nt = Th.nt;
    const int d = Mesh::Rd::d;
    const int dh = Mesh::RdHat::d;

    const long  nvk = nbvertexElement<Mesh>::value;
   if(verbosity>9)  cout << " nvk ="<< nvk << endl;
    if(pcc->N() != Th.nv) pcc->resize(Th.nv);
    KN<K> & cc=*pcc;
    long nbc =Th.nv;
    KN<long> forest(Th.nv,-1);// KrusKal algo
    for( int k=0; k< Th.nt;++k)
      for( int i1 =0; i1< nvk-1; ++i1 )// small opt because elem is connexe ..
       {
           long i2= (i1+1)%nvk;
           long s1 = Th(k,i1);
           long s2 = Th(k,i2);
           //  s0 and s1 are in the same cnx. comp.
           long r1=s1,r2=s2;
           // recherche des racines
           while( forest[r1]>=0) r1 = forest[r1];
           while( forest[r2]>=0) r2 = forest[r2];
           if(r1 !=r2)
           {
               nbc--;
               if(forest[r1]<forest[r2])
                   forest[r2]=r1; // lie  arbre r2 sur racine r1
               else if(forest[r2]<forest[r1])
                   forest[r1]=r2; // lie  arbre r2 sur racine r1
               else
                   forest[r1]=r2, forest[r2]--;  // lie  arbre r2 sur racine r1 et incremente al profondeur
           }
        }
    //  build cc
    cc=-1; //
    long nc =0;
    for(long s=0; s< nv; ++s)
    {
        long r=s;
        while( forest[r]>=0) r = forest[r];
        if( cc[r] <0) cc[r] = nc++;
        cc[s]=cc[r];
    }
    ffassert(nc==nbc);
    if( verbosity) cout << "  The number of  connexe componante (by vertex)  Mesh "<< pTh << " is " << nc << " / "
        " dim = " << d << " dim s " << dh <<  endl;
    return nc;
}

template<class Mesh,class K>
long connexecomponante(const Mesh * pTh,KN<K>* pcc,long flags)
{
    if(verbosity) cout << " ConnectedComponents closure flags "<< flags <<endl;
    long nbc =0;
    const Mesh & Th= *pTh;
    if(flags==1) //  data on element ..
    {
        KN<long> cv(pTh->nv);
        nbc = connexecomponantev(pTh,&cv);
        if( pcc->N() != Th.nv) pcc->resize(Th.nt);
        KN<K> & c = *pcc;
        const long  nvk = nbvertexElement<Mesh>::value;
        
        for( int k=0; k< Th.nt;++k)
          c[k]= cv[Th(k,0)];
    }
    else if (flags==2)
        nbc = connexecomponantev(pTh,pcc);
    else
        nbc = connexecomponantea(pTh,pcc);
    if(verbosity) cout << " nb. ConnectedComponents  "<< nbc << endl;
    return nbc;
}
template<class Mesh,class K>
class ConnectedComponents   :   public E_F0mps {
public:
 // typedef double K;
   typedef Mesh const * pmesh;
  typedef long  Result;
  Expression eTh,ecc;
  static basicAC_F0::name_and_type name_param[] ;
  static const int n_name_param =2;
  Expression nargs[n_name_param];

  long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}


  ConnectedComponents(const basicAC_F0 & args)
  {
      cout << "ConnectedComponents n_name_param" << n_name_param << endl;
    args.SetNameParam(n_name_param,name_param,nargs);
    eTh=to<pmesh>(args[0]);
    ecc=to<KN<K>*>(args[1]);
  }

  static ArrayOfaType  typeargs() {
    return  ArrayOfaType(atype<pmesh>(),atype<KN<K>*>(),false);}

  /// <<MeshCarre2_f>>

  static  E_F0 * f(const basicAC_F0 & args){
    return new ConnectedComponents(args);}

  AnyType operator()(Stack s) const {
    long flags=arg(0,s,0);
    flags = flags != 0;
    long flagsv=arg(1,s,0);
    if(flagsv) flags=2;
        pmesh pTh=  GetAny<pmesh>( (*eTh)(s));
         KN<K>* pcc = GetAny<KN<K>*>( (*ecc)(s));
    /// calls [[Carre]]

      return connexecomponante(pTh,pcc,flags);
  }

  operator aType () const { return atype<pmesh>();}
};

template<class Mesh,class Cmp=std::less<double> >
KN_<long> iminKP1(Stack stack,const Mesh * const & pTh,KN<double> * const & pu)
{
    const Cmp cmp;
    if(verbosity>9) cout << "iminKP1:  cmp(1.,2.) =" << cmp(1.,2.)<< endl;
    const int d= Mesh::RdHat::d;
    const int nbvK = d+1;
    ffassert(pu);
    const Mesh & Th=*pTh;
    KN<double> &u=*pu;
    ffassert(u.N()== Th.nv); // P1 ...
    long *pi = Add2StackOfPtr2FreeA(stack,new long[Th.nt]);
   if(verbosity>1) cout<< " i[min|max]KP1: nbvk ="<< nbvK << " nv " << Th.nv << " nt :" << Th.nt << " cmp: " << cmp(1.,2.) <<  endl;
    for (long k=0;k<Th.nt;++k)
    {
        long im=Th(k,0);
        for (int i=1;i<nbvK;++i)
        {
          long jm = Th(k,i);
 //           if( u[im] > u[jm] ) im =jm;
            if( cmp(u[jm], u[im]) ) im =jm;
         }
        pi[k]=im;
    }
    return KN_<long>(pi,Th.nt);
}
template<class Mesh >
KN_<long> imaxKP1(Stack stack,const Mesh * const & pTh,KN<double> * const & pu)
{ return iminKP1<Mesh,greater<double> >(stack,pTh,pu);}

template<class Mesh,class K>
basicAC_F0::name_and_type  ConnectedComponents<Mesh,K>::name_param[]= {
    {  "closure", &typeid(long) },
    {  "vertices", &typeid(long) }
};

static void inittt( ) {
    Global.Add("ConnectedComponents", "(",new OneOperatorCode<ConnectedComponents<Mesh,double> > );
    Global.Add("ConnectedComponents", "(",new OneOperatorCode<ConnectedComponents<Mesh,long> > );
    Global.Add("ConnectedComponents", "(",new OneOperatorCode<ConnectedComponents<Mesh3,double> > );
    Global.Add("ConnectedComponents", "(",new OneOperatorCode<ConnectedComponents<Mesh3,long> > );
    Global.Add("ConnectedComponents", "(",new OneOperatorCode<ConnectedComponents<MeshL,double> > );
    Global.Add("ConnectedComponents", "(",new OneOperatorCode<ConnectedComponents<MeshL,long> > );
    Global.Add("ConnectedComponents", "(",new OneOperatorCode<ConnectedComponents<MeshS,double> > );
    Global.Add("ConnectedComponents", "(",new OneOperatorCode<ConnectedComponents<MeshS,long> > );
    Global.Add("iminKP1", "(",new OneOperator2s_<KN_<long>,pmesh3,KN<double> * >(iminKP1));
    Global.Add("iminKP1", "(",new OneOperator2s_<KN_<long>,pmesh,KN<double> * >(iminKP1));
    Global.Add("iminKP1", "(",new OneOperator2s_<KN_<long>,pmeshL,KN<double> * >(iminKP1));
    Global.Add("iminKP1", "(",new OneOperator2s_<KN_<long>,pmeshS,KN<double> * >(iminKP1));
    Global.Add("imaxKP1", "(",new OneOperator2s_<KN_<long>,pmesh3,KN<double> * >(imaxKP1));
    Global.Add("imaxKP1", "(",new OneOperator2s_<KN_<long>,pmesh,KN<double> * >(imaxKP1));
    Global.Add("imaxKP1", "(",new OneOperator2s_<KN_<long>,pmeshL,KN<double> * >(imaxKP1));
    Global.Add("imaxKP1", "(",new OneOperator2s_<KN_<long>,pmeshS,KN<double> * >(imaxKP1));
/*
  Global.Add("connexecomponantev", "(",
             new OneOperator2<long,pmesh,KN<long>* >(connexecomponantev));
  Global.Add("connexecomponantev", "(",
               new OneOperator2<long,pmesh,KN<double>* >(connexecomponantev));
  Global.Add("connexecomponantev", "(",
               new OneOperator2<long,pmesh3,KN<long>* >(connexecomponantev));
  Global.Add("connexecomponantev", "(",
                 new OneOperator2<long,pmesh3,KN<double>* >(connexecomponantev));
    Global.Add("connexecomponantev", "(",
                 new OneOperator2<long,pmeshL,KN<long>* >(connexecomponantev));
    Global.Add("connexecomponantev", "(",
                   new OneOperator2<long,pmeshL,KN<double>* >(connexecomponantev));
    Global.Add("connexecomponantev", "(",
                 new OneOperator2<long,pmeshS,KN<long>* >(connexecomponantev));
    Global.Add("connexecomponantev", "(",
                   new OneOperator2<long,pmeshS,KN<double>* >(connexecomponantev));
 
    Global.Add("connexecomponantea", "(",
               new OneOperator2<long,pmesh,KN<long>* >(connexecomponantea));
    Global.Add("connexecomponantea", "(",
                 new OneOperator2<long,pmesh,KN<double>* >(connexecomponantea));

    Global.Add("connexecomponantea", "(",
                 new OneOperator2<long,pmesh3,KN<long>* >(connexecomponantea));
    Global.Add("connexecomponantea", "(",
                   new OneOperator2<long,pmesh3,KN<double>* >(connexecomponantea));
      Global.Add("connexecomponantea", "(",
                   new OneOperator2<long,pmeshL,KN<long>* >(connexecomponantea));
      Global.Add("connexecomponantea", "(",
                     new OneOperator2<long,pmeshL,KN<double>* >(connexecomponantea));
      Global.Add("connexecomponantea", "(",
                   new OneOperator2<long,pmeshS,KN<long>* >(connexecomponantea));
      Global.Add("connexecomponantea", "(",
                     new OneOperator2<long,pmeshS,KN<double>* >(connexecomponantea));
*/
}

LOADFUNC(inittt);
