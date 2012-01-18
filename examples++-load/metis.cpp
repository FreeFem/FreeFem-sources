//   for automatic  compilation with ff-c++
//ff-c++-LIBRARY-dep:   metis
//ff-c++-cpp-dep: 
//  
//   tools to read ppm file 
/*  use in freefem++ edp
  see :
  real[int,int] ff1("tt.pmm"); // read  image and set to an array. 
  real[int]  ff(ff1.nx*ff1.ny);
  ff=ff1; 
 */
#include <ff++.hpp>
#include <cmath>
typedef KNM<double> * pRnm;
typedef KN<double> * pRn;
typedef string * pstring;
extern "C" {
#include <metis.h>
}

#ifdef  METIS_VER_MAJOR 
//  METIS_PartMeshDual(&ne, &nn, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
extern "C" {
real_t libmetis__ComputeElementBalance(idx_t ne, idx_t nparts, idx_t *where); 
}
#else 
typedef idxtype idx_t ;
#endif
template<class Mesh,int NO>
KN<long> * partmetis(Stack s,KN<long> * const & part,Mesh * const & pTh,long const & lparts)
{
    ffassert(pTh);
    const Mesh & Th(*pTh);
    int nt=Th.nt,nv=Th.nv;
    int nve = Mesh::Rd::d+1;
    
    KN<idx_t> eptr(nt+1),elmnts(nve*nt), epart(nt), npart(nv);
    for(int k=0,i=0;k<nt;++k)
      {
	eptr[k]=i;
	for(int j=0;j<nve;j++)
	  elmnts[i++] = Th(k,j);
	eptr[k+1]=i;
      }
    int numflag=0;
    int nparts=lparts;
    int edgecut;
    int etype =nve-2; // triangle or tet .  change FH fevr 2010 
    idx_t ncommon = 1; 
#ifdef  METIS_VER_MAJOR 
    if(NO==0)
      METIS_PartMeshNodal(&nt, &nv, eptr, (idx_t *) elmnts,  0,0, &nparts, 0,0, &edgecut, (idx_t  *) epart, (idx_t  *) npart);
    else
      METIS_PartMeshDual(&nt, &nv, eptr, (idx_t *) elmnts , 0,0, &ncommon,  &nparts, 0,0, &edgecut, (idx_t  *) epart,(idx_t  *)  npart);
    if(verbosity)
      printf("  --metisOA: %d-way Edge-Cut: %7d, Balance: %5.2f Nodal=0/Dual %d\n", nparts, nve, libmetis__ComputeElementBalance(nt, nparts, epart),NO);
#else
    if(NO==0)
     METIS_PartMeshNodal(&nt, &nv, elmnts, &etype , &numflag, &nparts, &edgecut, epart, npart);
    else
     METIS_PartMeshDual(&nt, &nv, elmnts, &etype , &numflag, &nparts, &edgecut, epart, npart);
    if(verbosity)
      printf("  --metis: %d-way Edge-Cut: %7d, Balance: %5.2f Nodal=0/Dual %d\n", nparts, nve, ComputeElementBalance(nt, nparts, epart),NO);
#endif 
    part->resize(nt);
    *part=epart;
    return part;
}
KN<long> * partmetisd(Stack s,KN<long> * const & part,Mesh * const & pTh,long const & lparts)
{
    ffassert(pTh);
    const Mesh & Th(*pTh);
    int nt=Th.nt,nv=Th.nv;
    int nve = Mesh::Element::NbV;
    
    KN<idx_t> elmnts(nve*nt), epart(nt), npart(nv);
    for(int k=0,i=0;k<nt;++k)
	for(int j=0;j<nve;j++)
	    elmnts[i++] = Th(k,j);
    int numflag=0;
    int nparts=lparts;
    int edgecut;
    int etype =nve-2; // triangle
#ifdef  METIS_VER_MAJOR    
    printf("  %d-way Edge-Cut: %7d, Balance: %5.2f\n", nparts, nve, libmetis__ComputeElementBalance(nt, nparts, epart));
#else
    printf("  %d-way Edge-Cut: %7d, Balance: %5.2f\n", nparts, nve, ComputeElementBalance(nt, nparts, epart));
#endif
    part->resize(nt);
    *part=epart;
    return part;
}
class Init { public:
    Init();
};
// E_F_StackF0F0

LOADINIT(Init);
Init::Init(){
  if(verbosity && mpirank == 0) 
  cout << " lood: init metis  " << endl;
  Global.Add("metisnodal","(",new OneOperator3_<KN<long> *,KN<long> *,Mesh *,long , E_F_stackF0F0F0_<KN<long> *,KN<long> *,Mesh *,long> >(&partmetis<Mesh,0>));
  Global.Add("metisdual","(",new OneOperator3_<KN<long> *,KN<long> *,Mesh *,long , E_F_stackF0F0F0_<KN<long> *,KN<long> *,Mesh *,long> >(&partmetis<Mesh,1>));
    Global.Add("metisnodal","(",new OneOperator3_<KN<long> *,KN<long> *,Mesh3 *,long , E_F_stackF0F0F0_<KN<long> *,KN<long> *,Mesh3 *,long> >(&partmetis<Mesh3,0>));
    Global.Add("metisdual","(",new OneOperator3_<KN<long> *,KN<long> *,Mesh3 *,long , E_F_stackF0F0F0_<KN<long> *,KN<long> *,Mesh3 *,long> >(&partmetis<Mesh3,1>));
    
}
