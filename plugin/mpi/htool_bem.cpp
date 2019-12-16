//ff-c++-LIBRARY-dep: cxx11 [mkl|blas] mpi pthread htool bemtool boost
//ff-c++-cpp-dep:
// for def  M_PI under windows in <cmath>
#define _USE_MATH_DEFINES
#include <ff++.hpp>
#include <AFunction_ext.hpp>

#include <htool/lrmat/partialACA.hpp>
#include <htool/lrmat/fullACA.hpp>
#include <htool/lrmat/SVD.hpp>
#include <htool/types/matrix.hpp>
#include <htool/types/hmatrix.hpp>

// include the bemtool library .... path define in where library
//#include <bemtool/operator/block_op.hpp>  
#include <bemtool/tools.hpp>
#include <bemtool/fem/dof.hpp>
#include <bemtool/operator/operator.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
#include "PlotStream.hpp"

#include "common.hpp"

extern FILE *ThePlotStream;

using namespace std;
using namespace htool;
using namespace bemtool;

template <class MMesh>
class OrientNormal_Op : public E_F0mps {
 public:
  Expression eTh;
  static const int n_name_param = 1;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];

  bool arg(int i, Stack stack, bool a) const {
    return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
  }

 public:
  OrientNormal_Op(const basicAC_F0 &args, Expression tth) : eTh(tth) {
    args.SetNameParam(n_name_param, name_param, nargs);
  }
  AnyType operator( )(Stack stack) const;
};

inline double SolidAngle(const Fem2D::R3& p, const EdgeL& e){
  double M[2][2];
  for(int j=0; j<2; j++)
  for(int k=0;k<2;k++)
    M[j][k] = e[k][j] - p[j];
  return (1./2)*( M[0][0]*M[1][1]-M[1][0]*M[0][1] );
}

inline double SolidAngle(const Fem2D::R3& p, const TriangleS& e){
  double M[3][3];
  for(int j=0; j<3; j++)
  for(int k=0;k<3;k++)
    M[j][k] = e[k][j] - p[j];
  return (-1./6)*( M[0][0]*( M[1][1]*M[2][2]-M[2][1]*M[1][2] )
               - M[0][1]*( M[1][0]*M[2][2]-M[2][0]*M[1][2] )
               + M[0][2]*( M[1][0]*M[2][1]-M[2][0]*M[1][1] ) );
}

template <class MMesh>
basicAC_F0::name_and_type OrientNormal_Op<MMesh>::name_param[] = {

  {"unbounded", &typeid(bool)},
};

template <class Mesh>
void ComputeOrientation(const Mesh& Th, std::vector<bool>& orientation, bool unbounded) {
  typedef typename Mesh::RdHat RdHat;
  typedef typename Mesh::Element T;

  ////////////////// Compute connected component
  // tableau  avec les no. des elements
  // de chaque composante.
  // num[j][k] est le no du k ieme elt du
  // la composante no. j
  std::vector< std::vector<int> >  num;
  int nbc = 1;
  num.resize(nbc);
  int nb_visited = 0;
  std::queue<int> visit;
  int  nbelt = Th.nt;
  std::vector<bool> visited(nbelt,false);

  // Initialisation de l'algo
  int j0 = 0;
  visit.push(j0);
  visited[j0]=true;
  nb_visited++;
  num[nbc-1].push_back(j0);

  // Lancement de l'algo
  while(nb_visited<nbelt){
    // Reinitialisation dans le cas
    // de plusieurs composantes connexes
    if(visit.empty()){
      nbc++; num.resize(nbc);
      j0=0; while(visited[j0]){j0++;}
      visit.push(j0);
      if(!visited[j0]){
        visited[j0]=true;
        nb_visited++;
        num[nbc-1].push_back(j0);
      }
    }
    else{
      j0 = visit.front();
      visit.pop();
    }

    // Boucle sur les voisins de
    // l'element courant
    for(int k0=0; k0<RdHat::d+1; k0++){
      int k0a = k0;
      const int& j1 = Th.ElementAdj(j0, k0a);
      const T &K1(Th[j1]);
      if(!visited[j1]){
        nb_visited++;
        visit.push(j1);
        visited[j1]=true;
        num[nbc-1].push_back(j1);
      }
    }
  }

  ////////////////// Compute orientation
  // int  nbelt = Th.nt;
  int nt = Th.nt;
  bool ok = true;
  orientation.resize(nbelt,ok);
  std::fill(visited.begin(),visited.end(),false);
  std::vector<double> global_orientation(nbc);
  // std::vector<bool> visited(nbelt,false);

  //====================================
  // initialisation  de la recherche
  // d'un point extremal du maillage
  int  Iext = 0;
  RdHat bary = RdHat::diag(1./(RdHat::d+1));
  double Ext = Th[0](bary).norme2();
  //===============================//
  //   Breadth First Search sur    //
  //   chaque composante connexe   //
  //===============================//

  for(int I=0; I<nbc; I++){
    int nbe = num[I].size();
    int nb_visited = 0;
    std::queue<int> visit;

    int j0 = num[I][0];
    visit.push(j0);
    visited[j0]=true;

    while(nb_visited < nbe){

      j0 = visit.front();
      const T &K(Th[j0]);
      visit.pop();
      nb_visited++;

      if( K(bary).norme2() > Ext){
        Ext = K(bary).norme2();
        Iext = I;
      }

      for(int k0=0; k0<RdHat::d+1; k0++){
        int k0a = k0;
        const int& j1 = Th.ElementAdj(j0, k0a);
        const T &K1(Th[j1]);
        if(!visited[j1]){
          bool same = RdHat::d == 2 ? (K.EdgeOrientation(k0) != K1.EdgeOrientation(k0a)) : (k0 != k0a);
          if(same){orientation[j1]=orientation[j0];}
          else{orientation[j1]=!orientation[j0];}
          visited[j1]=true;
          visit.push(j1);
        }
      }
    }

    global_orientation[I] = 0.;
    const Fem2D::R3& p = Th.vertices[num[I][0]];
    for(int j=0; j<nbe; j++){
      j0 = num[I][j];
      const T &K1(Th[j0]);
      double r = SolidAngle(p,K1);
      if(!orientation[j0]){r = -r;}
      global_orientation[I] += r;
    }

    // for(int j=0; j<nt; j++) {
    //   if(global_orientation > 0)
    //     orientation[j] = !orientation[j];
    //   if (!unbounded)
    //     orientation[j] = !orientation[j];
    // }

    if(global_orientation[I] > 0){
      for(int j=0; j<nbe; j++) {
        j0 = num[I][j];
        orientation[j0] = !orientation[j0];
      }
    }
  }

  // //=====================================
  // // Calcul effectif des vecteurs normaux
  // normal.resize(nbelt);
  // for(int j=0; j<nbelt; j++){
  //   normal[j]=NormalTo(mesh[j]);
  //   normalize(normal[j]);
  //   if(orientation[j]){normal[j] = (-1.)*normal[j];}
  // }

  //=====================================
  // Si le domaine est borne la
  // composante exterieure du bord
  // doit etre orientee dans l'autre sens
  if(!unbounded){
    for(int j=0; j<num[Iext].size(); j++){
      int jj = num[Iext][j];
      orientation[j] = !orientation[j];
    }
  }
}

template <class MMesh>
AnyType OrientNormal_Op<MMesh>::operator( )(Stack stack) const {
  typedef typename MMesh::Element T;
  typedef typename MMesh::BorderElement B;
  typedef typename MMesh::Vertex V;

  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  MMesh *pTh = GetAny< MMesh * >((*eTh)(stack));
  MMesh &Th = *pTh;
  ffassert(pTh);

  bool unbounded(arg(0, stack, false));

  if (verbosity > 5) cout << "Orienting surface normals ..." << endl;

  int nv = Th.nv, nt = Th.nt, nbe = Th.nbe;

  V *v = new V[nv];
  T *t = new T[nt];
  T *tt = t;
  B *b = new B[nbe];
  B *bb = b;
  double mes = 0, mesb = 0;

  if (verbosity > 5)
    cout << "copy the original mesh ... nv= " << nv << " nt= " << nt << " nbe= " << nbe << endl;

  for (int i = 0; i < nv; i++) {
    const V &K(Th.vertices[i]);
    v[i].x = K.x;
    v[i].y = K.y;
    v[i].z = K.z;
    v[i].lab = K.lab;
  }

  std::vector<bool> orientation;
  ComputeOrientation(Th, orientation, unbounded);
  //=====================================
  // Add elements, with correct orientation

  for (int i = 0; i < nt; i++) {
    const T &K(Th.elements[i]);
    int iv[T::nea];
    int lab = K.lab;

    for (int jj = 0; jj < T::nea; jj++) {
      iv[jj] = Th.operator( )(K[jj]);
      assert(iv[jj] >= 0 && iv[jj] < nv);
    }
    if (!orientation[i])
        swap(iv[0], iv[1]);
    (tt)->set(v, iv, lab);
    mes += tt++->mesure();
  }

  for (int i = 0; i < nbe; i++) {
    const B &K(Th.be(i));
    int iv[B::nea];
    int lab = K.lab;
    for (int jj = 0; jj < B::nea; jj++) {
      iv[jj] = Th.operator( )(K[jj]);
      assert(iv[jj] >= 0 && iv[jj] < nv);
    }
    (bb)->set(v, iv, lab);
    mesb += bb++->mesure( );
  }

  MMesh *Th_t = new MMesh(nv, nt, nbe, v, t, b);
  Th_t->BuildGTree( );
  //Th_t->BuildMeshL( );
  *mp = mps;
  Add2StackOfPtr2FreeRC(stack, Th_t);
  return Th_t;
}

template <class MMesh>
class OrientNormal : public OneOperator {
 public:
  OrientNormal( ) : OneOperator(atype< const MMesh* >( ), atype< const MMesh* >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new OrientNormal_Op<MMesh>(args, t[0]->CastTo(args[0]));
  }
};

template<class K>
class MyMatrix: public IMatrix<K>{
	const MeshS & ThU; // line
	const MeshS & ThV; // colunm

public:
	MyMatrix(const FESpaceS * Uh , const FESpaceS * Vh ):IMatrix<K>(Uh->Th.nv,Vh->Th.nv),ThU(Uh->Th), ThV(Vh->Th) {}

	K get_coef(const int& i, const int& j)const {return 1./(0.01+(ThU.vertices[i]-ThV.vertices[j]).norme2());}

};

template<class K>
class HMatrixVirt {
public:
		virtual const std::map<std::string, std::string>& get_infos() const = 0;
		virtual void mvprod_global(const K* const in, K* const out,const int& mu=1) const = 0;
		virtual int nb_rows() const = 0;
		virtual int nb_cols() const = 0;
		virtual void cluster_to_target_permutation(const K* const in, K* const out) const = 0;
		virtual const MPI_Comm& get_comm() const = 0;
		virtual int get_rankworld() const = 0;
		virtual int get_sizeworld() const = 0;
		virtual const std::vector<SubMatrix<K>*>& get_MyNearFieldMats() const = 0;
		virtual const LowRankMatrix<K>& get_MyFarFieldMats(int i) const = 0;
		virtual int get_MyFarFieldMats_size() const = 0;
		virtual const std::vector<SubMatrix<K>*>& get_MyStrictlyDiagNearFieldMats() const = 0;
		virtual Matrix<K> to_dense_perm() const = 0;

		virtual ~HMatrixVirt() {};
};

template<template<class> class LR, class K>
class HMatrixImpl : public HMatrixVirt<K> {
private:
		HMatrix<LR,K> H;
public:
		HMatrixImpl(IMatrix<K>& I, const std::vector<htool::R3>& xt, const int& reqrank=-1, MPI_Comm comm=MPI_COMM_WORLD) : H(I,xt,reqrank,comm){}
		HMatrixImpl(IMatrix<K>& I, const std::vector<htool::R3>& xt, const std::vector<htool::R3>& xs, const int& reqrank=-1, MPI_Comm comm=MPI_COMM_WORLD) : H(I,xt,xs,reqrank,comm){}
		const std::map<std::string, std::string>& get_infos() const {return H.get_infos();}
		void mvprod_global(const K* const in, K* const out,const int& mu=1) const {return H.mvprod_global(in,out,mu);}
		int nb_rows() const { return H.nb_rows();}
		int nb_cols() const { return H.nb_cols();}
		void cluster_to_target_permutation(const K* const in, K* const out) const {return H.cluster_to_target_permutation(in,out);}
		const MPI_Comm& get_comm() const {return H.get_comm();}
		int get_rankworld() const {return H.get_rankworld();}
		int get_sizeworld() const {return H.get_sizeworld();}
		const std::vector<SubMatrix<K>*>& get_MyNearFieldMats() const {return H.get_MyNearFieldMats();}
		const LowRankMatrix<K>& get_MyFarFieldMats(int i) const {return *(H.get_MyFarFieldMats()[i]);}
		int get_MyFarFieldMats_size() const {return H.get_MyFarFieldMats().size();}
		const std::vector<SubMatrix<K>*>& get_MyStrictlyDiagNearFieldMats() const {return H.get_MyStrictlyDiagNearFieldMats();}
		Matrix<K> to_dense_perm() const {return H.to_dense_perm();}
};

template<class v_fes1, class v_fes2, class K>
class assembleHMatrix : public OneOperator { public:

	class Op : public E_F0info {
	public:
		Expression a,b,c;

		static const int n_name_param = 10;
		static basicAC_F0::name_and_type name_param[] ;
		Expression nargs[n_name_param];
		bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
		long argl(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
		string* args(int i,Stack stack,string* a) const{ return nargs[i] ? GetAny<string*>( (*nargs[i])(stack) ): a;}
		double arg(int i,Stack stack,double a) const{ return nargs[i] ? GetAny<double>( (*nargs[i])(stack) ): a;}
		K argk(int i,Stack stack,K a) const{ return nargs[i] ? GetAny<K>( (*nargs[i])(stack) ): a;}
		KN_<long> arg(int i,Stack stack,KN_<long> a ) const{ return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
		pcommworld arg(int i,Stack stack,pcommworld a ) const{ return nargs[i] ? GetAny<pcommworld>( (*nargs[i])(stack) ): a;}
	public:
		Op(const basicAC_F0 &  args,Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc) {
			args.SetNameParam(n_name_param,name_param,nargs); }
		};

		assembleHMatrix() : OneOperator(atype<const typename assembleHMatrix<v_fes1,v_fes2,K>::Op *>(),
		atype<v_fes1 **>(),
		atype<v_fes2 **>(),
		atype<K*>()) {}

		E_F0 * code(const basicAC_F0 & args) const
		{
			return  new Op(args,t[0]->CastTo(args[0]),
			t[1]->CastTo(args[1]),
			t[2]->CastTo(args[2]));
		}
};

template<class v_fes1, class v_fes2, class K>
basicAC_F0::name_and_type  assembleHMatrix<v_fes1,v_fes2,K>::Op::name_param[]= {
		{  "eps", &typeid(double)},
		{  "eta", &typeid(double)},
		{  "minclustersize", &typeid(long)},
		{  "maxblocksize", &typeid(long)},
		{  "mintargetdepth", &typeid(long)},
		{  "minsourcedepth", &typeid(long)},
		{  "commworld", &typeid(pcommworld)},
		{  "alpha", &typeid(double)},
		{  "compressor", &typeid(string*)},
		{  "combinedcoef", &typeid(K)}
};

template<class ffmesh, class bemtoolmesh>
void Mesh2Bemtool(const ffmesh &Th, Geometry &node, bemtoolmesh &mesh ) {
   
    typedef typename ffmesh::RdHat RdHat;
    typedef typename ffmesh::Element E;
    const int dHat =  RdHat::d;
    
    // create the geometry;
    
    bemtool::R3 p;
    for(int iv=0 ; iv<Th.nv ; iv++){
        p[0]=Th.vertices[iv].x;p[1]=Th.vertices[iv].y;p[2]=Th.vertices[iv].z;
        node.setnodes(p);
    }
   
    node.initEltData();
   
    if(verbosity>10) std::cout << "Creating mesh domain (nodes)" << std::endl;
  
    mesh.set_elt(node);
    bemtool::array<dHat+1,int> I;
    if(verbosity>10) std::cout << "End creating mesh domain mesh" << std::endl;

    if(verbosity>10) std::cout << "Creating geometry domain (elements)" << std::endl;
    for(int it=0; it<Th.nt; it++){
      const E &K(Th[it]);
      for(int j=0;j<dHat+1;j++)
        I[j]=Th.operator () (K[j]);
      mesh.setOneElt(node,I);
    }

    //mesh = unbounded;
    //Orienting(mesh);
    Normal<dHat> N(mesh);
    for(int it=0; it<Th.nt; it++){
      const E &K(Th[it]);
      Fem2D::R3 nn = K.NormalSUnitaire();
      bemtool::R3 mm; mm[0]=nn.x; mm[1]=nn.y; mm[2]=nn.z;
      N.set(it, mm);
    }
    mesh.Orienting(N);

    if(verbosity>10) std::cout << "end creating geometry domain" << std::endl;
}

template<class ffmesh>
void Mesh2Bemtool(const ffmesh &Th, Geometry &node) {
    if(verbosity>10) std::cout << "Creating mesh output" << std::endl;
    bemtool::R3 p;
    Fem2D::R3 pp;
    for(int iv=0 ; iv<Th.nv ; iv++){
      pp = Th.vertices[iv];
      p[0]=pp.x; p[1]=pp.y; p[2]=pp.z;
      node.setnodes(p);
    }
}

enum TypeEnum { tBIOp, tBIOpwithmass, tBIOpcombined, tPotential };

template<class K, class v_fes1, class v_fes2, EquationEnum Eq, TypeEnum type, BIOpKernelEnum Ker, BIOpKernelEnum Ker2, PotKernelEnum Pot>
AnyType SetHMatrix(Stack stack,Expression emat,Expression einter,int init)
{
	using namespace Fem2D;

	HMatrixVirt<K>** Hmat =GetAny<HMatrixVirt<K>** >((*emat)(stack));
	const typename assembleHMatrix<v_fes1,v_fes2,K>::Op * mi(dynamic_cast<const typename assembleHMatrix<v_fes1,v_fes2,K>::Op *>(einter));

	double epsilon=mi->arg(0,stack,1e-2);//htool::Parametres::epsilon);
	double eta=mi->arg(1,stack,10.);//htool::Parametres::eta);
	int minclustersize=mi->argl(2,stack,10);//htool::Parametres::minclustersize);
	int maxblocksize=mi->argl(3,stack,1000000);//htool::Parametres::maxblocksize);
	int mintargetdepth=mi->argl(4,stack,1);//htool::Parametres::mintargetdepth);
	int minsourcedepth=mi->argl(5,stack,1);//htool::Parametres::minsourcedepth);
	pcommworld pcomm=mi->arg(6,stack,nullptr);
	double alpha=mi->arg(7,stack,0.5);
	string* compressor=mi->args(8,stack,0);
	K combinedcoef=mi->argk(9,stack,0.5);

	SetMaxBlockSize(maxblocksize);
	SetMinClusterSize(minclustersize);
	SetEpsilon(epsilon);
	SetEta(eta);
	SetMinTargetDepth(mintargetdepth);
	SetMinSourceDepth(minsourcedepth);

	MPI_Comm comm = pcomm ? *(MPI_Comm*)pcomm : MPI_COMM_WORLD;

	ffassert(einter);

	typedef typename v_fes1::FESpace FESpace1;
	typedef typename v_fes2::FESpace FESpace2;
	typedef typename FESpace1::Mesh Mesh1;
	typedef typename FESpace2::Mesh Mesh2;

	typedef typename std::conditional<Mesh1::RdHat::d==1,Mesh1D,Mesh2D>::type Mesh1Bemtool;
	typedef typename std::conditional<Mesh1::RdHat::d==1,P1_1D,P1_2D>::type P1;

	v_fes1** pUh = GetAny<v_fes1**>((* mi->a)(stack));
	FESpace1 * Uh = **pUh;
	int NUh =Uh->N;

	v_fes2** pVh = GetAny<v_fes2**>((* mi->b)(stack));
	FESpace2 * Vh = **pVh;
	int NVh =Vh->N;

	K * coef = GetAny< K * >((* mi->c)(stack));
	double kappa = coef->real();

	ffassert(Vh);
	ffassert(Uh);

	int n=Uh->NbOfDF;
	int m=Vh->NbOfDF;

	const Mesh1 & ThU =Uh->Th; // line
	const Mesh2 & ThV =Vh->Th; // colunm

	bool samemesh = (void*)&Uh->Th == (void*)&Vh->Th;  // same Fem2D::Mesh

	bool isPot = (type == tPotential);
	if (!isPot) ffassert (samemesh);

	if (init)
	  delete *Hmat;

	// call interface MeshS2Bemtool
	Geometry node; Mesh1Bemtool mesh;
	Mesh2Bemtool(ThU, node, mesh);

	if(verbosity>10) std::cout << "Creating dof" << std::endl;
	Dof<P1> dof(mesh,true);
	// now the list of dof is known -> can acces to global num triangle and the local num vertice assiciated 

	vector<htool::R3> p1(n);
	vector<htool::R3> p2(m);

	Fem2D::R3 pp;
	for (int i=0; i<n; i++) {
	  pp = ThU.vertices[i];
	  p1[i] = {pp.x, pp.y, pp.z};
	}

	if(!samemesh) {
	  for (int i=0; i<m; i++) {
	    pp = ThV.vertices[i];
	    p2[i] = {pp.x, pp.y, pp.z};
	  }
	}
	else
	  p2=p1;

	if (!isPot) {
	  // BIO_Generator
	  IMatrix<K>* generator;

	  if (type == tBIOp)
	    generator = new BIO_Generator<BIOpKernel<Eq,Ker,Mesh1::RdHat::d+1,P1,P1>,P1>(dof,kappa);
	  else if (type == tBIOpwithmass)
	    generator = new BIO_Generator_w_mass<BIOpKernel<Eq,Ker,Mesh1::RdHat::d+1,P1,P1>,P1>(dof,kappa,alpha);
	  else if (type == tBIOpcombined)
	    generator = new Combined_BIO_Generator<BIOpKernel<Eq,Ker,Mesh1::RdHat::d+1,P1,P1>,BIOpKernel<Eq,Ker2,Mesh1::RdHat::d+1,P1,P1>,P1>(dof,kappa,combinedcoef,alpha);

	  if (!compressor || *compressor == "partialACA")
	    *Hmat = new HMatrixImpl<partialACA,K>(*generator,p1,-1,comm);
	  else if (*compressor == "fullACA")
	    *Hmat = new HMatrixImpl<fullACA,K>(*generator,p1,-1,comm);
	  else if (*compressor == "SVD")
	    *Hmat = new HMatrixImpl<SVD,K>(*generator,p1,-1,comm);
	  else {
	    cerr << "Error: unknown htool compressor \""+*compressor+"\"" << endl;
	    ffassert(0);
	  }
	  delete generator;
	}
	else {
	  Geometry node_output;
	  Mesh2Bemtool(ThV,node_output);

	  Potential<PotKernel<Eq,Pot,Mesh1::RdHat::d+1,P1>> P(mesh,kappa);
	  POT_Generator<PotKernel<Eq,Pot,Mesh1::RdHat::d+1,P1>,P1> generator(P,dof,node_output);

	  if (!compressor || *compressor == "partialACA")
	    *Hmat = new HMatrixImpl<partialACA,K>(generator,p2,p1,-1,comm);
	  else if (*compressor == "fullACA")
	    *Hmat = new HMatrixImpl<fullACA,K>(generator,p2,p1,-1,comm);
	  else if (*compressor == "SVD")
	    *Hmat = new HMatrixImpl<SVD,K>(generator,p2,p1,-1,comm);
	  else {
	    cerr << "Error: unknown htool compressor \""+*compressor+"\"" << endl;
	    ffassert(0);
	  }
	}

	return Hmat;
}

template<class K, class v_fes1, class v_fes2, EquationEnum Eq,  TypeEnum type, BIOpKernelEnum Ker, BIOpKernelEnum Ker2, PotKernelEnum Pot, int init>
AnyType SetHMatrix(Stack stack,Expression emat,Expression einter)
{ return SetHMatrix<K,v_fes1,v_fes2,Eq,type,Ker,Ker2,Pot>(stack,emat,einter,init);}

template<class K>
AnyType ToDense(Stack stack,Expression emat,Expression einter,int init)
{
	ffassert(einter);
	HMatrixVirt<K>** Hmat =GetAny<HMatrixVirt<K>** >((*einter)(stack));
	ffassert(Hmat && *Hmat);
	HMatrixVirt<K>& H = **Hmat;
	Matrix<K> mdense = H.to_dense_perm();
	const std::vector<K>& vdense = mdense.get_mat();

	KNM<K>* M =GetAny<KNM<K>*>((*emat)(stack));

	for (int i=0; i< mdense.nb_rows(); i++)
		for (int j=0; j< mdense.nb_cols(); j++)
			(*M)(i,j) = mdense(i,j);

	return M;
}

template<class K, int init>
AnyType ToDense(Stack stack,Expression emat,Expression einter)
{ return ToDense<K>(stack,emat,einter,init);}

template<class V, class K>
class Prod {
public:
	const HMatrixVirt<K>* h;
	const V u;
	Prod(HMatrixVirt<K>** v, V w) : h(*v), u(w) {}

	void prod(V x) const {h->mvprod_global(*(this->u), *x);};

	static V mv(V Ax, Prod<V, K> A) {
		*Ax = K();
		A.prod(Ax);
		return Ax;
	}
	static V init(V Ax, Prod<V, K> A) {
		Ax->init(A.u->n);
		return mv(Ax, A);
	}

};

template<class K>
std::map<std::string, std::string>* get_infos(HMatrixVirt<K>** const& H) {
	return new std::map<std::string, std::string>((*H)->get_infos());
}

string* get_info(std::map<std::string, std::string>* const& infos, string* const& s){
	return new string((*infos)[*s]);
}

ostream & operator << (ostream &out, const std::map<std::string, std::string> & infos)
{
	for (std::map<std::string,std::string>::const_iterator it = infos.begin() ; it != infos.end() ; ++it){
		out<<it->first<<"\t"<<it->second<<std::endl;
	}
	out << std::endl;
	return out;
}

template<class A>
struct PrintPinfos: public binary_function<ostream*,A,ostream*> {
	static ostream* f(ostream* const  & a,const A & b)  {  *a << *b;
		return a;
	}
};

template<class K>
class plotHMatrix : public OneOperator {
public:

	class Op : public E_F0info {
	public:
		Expression a;

		static const int n_name_param = 2;
		static basicAC_F0::name_and_type name_param[] ;
		Expression nargs[n_name_param];
		bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
		long argl(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}

	public:
		Op(const basicAC_F0 &  args,Expression aa) : a(aa) {
			args.SetNameParam(n_name_param,name_param,nargs);
		}

		AnyType operator()(Stack stack) const{

			bool wait = arg(0,stack,false);
			long dim = argl(1,stack,2);

			HMatrixVirt<K>** H =GetAny<HMatrixVirt<K>** >((*a)(stack));

			PlotStream theplot(ThePlotStream);

			if (mpirank == 0 && ThePlotStream) {
				theplot.SendNewPlot();
				theplot << 3L;
				theplot <= wait;
				theplot << 17L;
				theplot <= dim;
				theplot.SendEndArgPlot();
				theplot.SendMeshes();
				theplot << 0L;
				theplot.SendPlots();
				theplot << 1L;
				theplot << 31L;
			}

			if (!H || !(*H)) {
				if (mpirank == 0&& ThePlotStream) {
					theplot << 0;
					theplot << 0;
					theplot << 0L;
					theplot << 0L;
				}
			}
			else {
				const std::vector<SubMatrix<K>*>& dmats = (*H)->get_MyNearFieldMats();

				int nbdense = dmats.size();
				int nblr = (*H)->get_MyFarFieldMats_size();

				int sizeworld = (*H)->get_sizeworld();
				int rankworld = (*H)->get_rankworld();

				int nbdenseworld[sizeworld];
				int nblrworld[sizeworld];
				MPI_Allgather(&nbdense, 1, MPI_INT, nbdenseworld, 1, MPI_INT, (*H)->get_comm());
				MPI_Allgather(&nblr, 1, MPI_INT, nblrworld, 1, MPI_INT, (*H)->get_comm());
				int nbdenseg = 0;
				int nblrg = 0;
				for (int i=0; i<sizeworld; i++) {
					nbdenseg += nbdenseworld[i];
					nblrg += nblrworld[i];
				}

				int* buf = new int[4*(mpirank==0?nbdenseg:nbdense) + 5*(mpirank==0?nblrg:nblr)];

				for (int i=0;i<nbdense;i++) {
					const SubMatrix<K>& l = *(dmats[i]);
					buf[4*i] = l.get_offset_i();
					buf[4*i+1] = l.get_offset_j();
					buf[4*i+2] = l.nb_rows();
					buf[4*i+3] = l.nb_cols();
				}

				int displs[sizeworld];
				int recvcounts[sizeworld];
				displs[0] = 0;

				for (int i=0; i<sizeworld; i++) {
					recvcounts[i] = 4*nbdenseworld[i];
					if (i > 0)	displs[i] = displs[i-1] + recvcounts[i-1];
				}
				MPI_Gatherv(rankworld==0?MPI_IN_PLACE:buf, recvcounts[rankworld], MPI_INT, buf, recvcounts, displs, MPI_INT, 0, (*H)->get_comm());

				int* buflr = buf + 4*(mpirank==0?nbdenseg:nbdense);
				double* bufcomp = new double[mpirank==0?nblrg:nblr];

				for (int i=0;i<nblr;i++) {
					const LowRankMatrix<K>& l = (*H)->get_MyFarFieldMats(i);
					buflr[5*i] = l.get_offset_i();
					buflr[5*i+1] = l.get_offset_j();
					buflr[5*i+2] = l.nb_rows();
					buflr[5*i+3] = l.nb_cols();
					buflr[5*i+4] = l.rank_of();
					bufcomp[i] = l.compression();
				}

				for (int i=0; i<sizeworld; i++) {
					recvcounts[i] = 5*nblrworld[i];
					if (i > 0)	displs[i] = displs[i-1] + recvcounts[i-1];
				}

				MPI_Gatherv(rankworld==0?MPI_IN_PLACE:buflr, recvcounts[rankworld], MPI_INT, buflr, recvcounts, displs, MPI_INT, 0, (*H)->get_comm());

				for (int i=0; i<sizeworld; i++) {
					recvcounts[i] = nblrworld[i];
					if (i > 0)	displs[i] = displs[i-1] + recvcounts[i-1];
				}

				MPI_Gatherv(rankworld==0?MPI_IN_PLACE:bufcomp, recvcounts[rankworld], MPI_DOUBLE, bufcomp, recvcounts, displs, MPI_DOUBLE, 0, (*H)->get_comm());

				if (mpirank == 0 && ThePlotStream ) {

					int si = (*H)->nb_rows();
					int sj = (*H)->nb_cols();

					theplot << si;
					theplot << sj;
					theplot << (long)nbdenseg;
					theplot << (long)nblrg;

					for (int i=0;i<nbdenseg;i++) {
						theplot << buf[4*i];
						theplot << buf[4*i+1];
						theplot << buf[4*i+2];
						theplot << buf[4*i+3];
					}

					for (int i=0;i<nblrg;i++) {
						theplot << buflr[5*i];
						theplot << buflr[5*i+1];
						theplot << buflr[5*i+2];
						theplot << buflr[5*i+3];
						theplot << buflr[5*i+4];
						theplot << bufcomp[i];
					}

					theplot.SendEndPlot();

				}
				delete [] buf;
				delete [] bufcomp;

			}

			return 0L;
		}
	};

	plotHMatrix() : OneOperator(atype<long>(),atype<HMatrixVirt<K> **>()) {}

	E_F0 * code(const basicAC_F0 & args) const
	{
		return  new Op(args,t[0]->CastTo(args[0]));
	}
};

template<class K>
basicAC_F0::name_and_type  plotHMatrix<K>::Op::name_param[]= {
	{  "wait", &typeid(bool)},
	{  "dim", &typeid(long)}
};

template<class T, class U, class K, char trans>
class HMatrixInv {
    public:
        const T t;
        const U u;

        struct HMatVirt: CGMatVirt<int,K> {
            const T tt;

            HMatVirt(T ttt) : tt(ttt), CGMatVirt<int,K>((*ttt)->nb_rows()) {}
            K*  addmatmul(K* x,K* Ax) const { (*tt)->mvprod_global(x, Ax); return Ax;}
        };

        struct HMatVirtPrec: CGMatVirt<int,K> {
            const T tt;
            std::vector<K> invdiag;

            HMatVirtPrec(T ttt) : tt(ttt), CGMatVirt<int,K>((*ttt)->nb_rows()), invdiag((*ttt)->nb_rows(),0) {
              std::vector<SubMatrix<K>*> diagblocks = (*tt)->get_MyStrictlyDiagNearFieldMats();
              std::vector<K> tmp((*ttt)->nb_rows(),0);
              for (int j=0;j<diagblocks.size();j++){
                SubMatrix<K>& submat = *(diagblocks[j]);
                int local_nr = submat.nb_rows();
                int local_nc = submat.nb_cols();
                int offset_i = submat.get_offset_i();
                int offset_j = submat.get_offset_j();
                for (int i=offset_i;i<offset_i+std::min(local_nr,local_nc);i++){
                  tmp[i] = 1./submat(i-offset_i,i-offset_i);
                }
              }
            (*tt)->cluster_to_target_permutation(tmp.data(),invdiag.data());
            MPI_Allreduce(MPI_IN_PLACE, &(invdiag[0]), (*ttt)->nb_rows(), wrapper_mpi<K>::mpi_type(), MPI_SUM, (*tt)->get_comm());
            }

            K*  addmatmul(K* x,K* Ax) const {
              for (int i=0; i<(*tt)->nb_rows(); i++)
                Ax[i] = invdiag[i] * x[i];
              return Ax;
            }
        };

        HMatrixInv(T v, U w) : t(v), u(w) {}

        void solve(U out) const {
            HMatVirt A(t);
            HMatVirtPrec P(t);
            bool res=fgmres(A,P,1,(K*)*u,(K*)*out,1.e-6,2000,200,(mpirank==0)*verbosity);
        }

        static U inv(U Ax, HMatrixInv<T, U, K, trans> A) {
            A.solve(Ax);
            return Ax;
        }
        static U init(U Ax, HMatrixInv<T, U, K, trans> A) {
            Ax->init(A.u->n);
            return inv(Ax, A);
        }
};

template<class K>
void addHmat() {
	Dcl_Type<HMatrixVirt<K>**>(Initialize<HMatrixVirt<K>*>, Delete<HMatrixVirt<K>*>);
	Dcl_TypeandPtr<HMatrixVirt<K>*>(0,0,::InitializePtr<HMatrixVirt<K>*>,::DeletePtr<HMatrixVirt<K>*>);
	//atype<HMatrix<LR ,K>**>()->Add("(","",new OneOperator2_<string*, HMatrix<LR ,K>**, string*>(get_infos<LR,K>));
	Add<HMatrixVirt<K>**>("infos",".",new OneOperator1_<std::map<std::string, std::string>*, HMatrixVirt<K>**>(get_infos));

	Dcl_Type<Prod<KN<K>*, K>>();
	TheOperators->Add("*", new OneOperator2<Prod<KN<K>*, K>, HMatrixVirt<K>**, KN<K>*>(Build));
	TheOperators->Add("=", new OneOperator2<KN<K>*, KN<K>*, Prod<KN<K>*, K>>(Prod<KN<K>*, K>::mv));
	TheOperators->Add("<-", new OneOperator2<KN<K>*, KN<K>*, Prod<KN<K>*, K>>(Prod<KN<K>*, K>::init));

	addInv<HMatrixVirt<K>*, HMatrixInv, KN<K>, K>();

	Global.Add("display","(",new plotHMatrix<K>);

	// to dense:
	TheOperators->Add("=",
	new OneOperator2_<KNM<K>*, KNM<K>*, HMatrixVirt<K>**,E_F_StackF0F0>(ToDense<K, 1>));
	TheOperators->Add("<-",
	new OneOperator2_<KNM<K>*, KNM<K>*, HMatrixVirt<K>**,E_F_StackF0F0>(ToDense<K, 0>));
}

template<class K, class v_fes1, class v_fes2, EquationEnum Eq, TypeEnum type, BIOpKernelEnum Ker, BIOpKernelEnum Ker2 = SL_OP>
void addBIOp(const char* namec) {
	Dcl_Type<const typename assembleHMatrix<v_fes1,v_fes2,K>::Op *>();
	Add<const typename assembleHMatrix<v_fes1,v_fes2,K>::Op *>("<-","(", new assembleHMatrix<v_fes1,v_fes2,K>);

	TheOperators->Add("=",
	new OneOperator2_<HMatrixVirt<K>**,HMatrixVirt<K>**,const typename assembleHMatrix<v_fes1,v_fes2,K>::Op*,E_F_StackF0F0>(SetHMatrix<K,v_fes1,v_fes2,Eq,type,Ker,Ker2,SL_POT, 1>));
	TheOperators->Add("<-",
	new OneOperator2_<HMatrixVirt<K>**,HMatrixVirt<K>**,const typename assembleHMatrix<v_fes1,v_fes2,K>::Op*,E_F_StackF0F0>(SetHMatrix<K,v_fes1,v_fes2,Eq,type,Ker,Ker2,SL_POT, 0>));

	Global.Add(namec,"(",new assembleHMatrix<v_fes1,v_fes2,K>);
}

template<class K, class v_fes1, class v_fes2, EquationEnum Eq, PotKernelEnum Pot>
void addPotential(const char* namec) {
	Dcl_Type<const typename assembleHMatrix<v_fes1,v_fes2,K>::Op *>();
	Add<const typename assembleHMatrix<v_fes1,v_fes2,K>::Op *>("<-","(", new assembleHMatrix<v_fes1,v_fes2,K>);

	TheOperators->Add("=",
	new OneOperator2_<HMatrixVirt<K>**,HMatrixVirt<K>**,const typename assembleHMatrix<v_fes1,v_fes2,K>::Op*,E_F_StackF0F0>(SetHMatrix<K,v_fes1,v_fes2,Eq,tPotential,SL_OP,SL_OP,Pot, 1>));
	TheOperators->Add("<-",
	new OneOperator2_<HMatrixVirt<K>**,HMatrixVirt<K>**,const typename assembleHMatrix<v_fes1,v_fes2,K>::Op*,E_F_StackF0F0>(SetHMatrix<K,v_fes1,v_fes2,Eq,tPotential,SL_OP,SL_OP,Pot, 0>));

	Global.Add(namec,"(",new assembleHMatrix<v_fes1,v_fes2,K>);
}

static void Init_Bem() {
	Global.Add("OrientNormal", "(", new OrientNormal<MeshS>);
	Global.Add("OrientNormal", "(", new OrientNormal<MeshL>);

	Dcl_Type<std::map<std::string, std::string>*>( );
	TheOperators->Add("<<",new OneBinaryOperator<PrintPinfos<std::map<std::string, std::string>*>>);
	Add<std::map<std::string, std::string>*>("[","",new OneOperator2_<string*, std::map<std::string, std::string>*, string*>(get_info));

	addHmat<std::complex<double>>();
	//add<partialACA,double>("assemble");

	// 3D BIOp
	addBIOp<std::complex<double>,v_fesS,v_fesS,HE,tBIOp,SL_OP>("assemblecomplexHESL");
	addBIOp<std::complex<double>,v_fesS,v_fesS,HE,tBIOp,DL_OP>("assemblecomplexHEDL");
	addBIOp<std::complex<double>,v_fesS,v_fesS,HE,tBIOp,HS_OP>("assemblecomplexHEHS");
	addBIOp<std::complex<double>,v_fesS,v_fesS,HE,tBIOpwithmass,DL_OP>("assemblecomplexHEDLwmass");
	addBIOp<std::complex<double>,v_fesS,v_fesS,HE,tBIOpwithmass,TDL_OP>("assemblecomplexHETDLwmass");

	addBIOp<std::complex<double>,v_fesS,v_fesS,HE,tBIOpcombined,SL_OP,DL_OP>("assemblecomplexHEcombinedSLDL");
	addBIOp<std::complex<double>,v_fesS,v_fesS,HE,tBIOpcombined,SL_OP,TDL_OP>("assemblecomplexHEcombinedSLTDL");
	addBIOp<std::complex<double>,v_fesS,v_fesS,HE,tBIOpcombined,HS_OP,DL_OP>("assemblecomplexHEcombinedHSDL");
	addBIOp<std::complex<double>,v_fesS,v_fesS,HE,tBIOpcombined,HS_OP,TDL_OP>("assemblecomplexHEcombinedHSTDL");

	// 3D Potential
	addPotential<std::complex<double>,v_fesS,v_fesS,HE,SL_POT>("assemblecomplexHESLPot");
	addPotential<std::complex<double>,v_fesS,v_fesS,HE,DL_POT>("assemblecomplexHEDLPot");

	// 2D BIOp
	addBIOp<std::complex<double>,v_fesL,v_fesL,HE,tBIOp,SL_OP>("assemblecomplexHESL");
	addBIOp<std::complex<double>,v_fesL,v_fesL,HE,tBIOp,DL_OP>("assemblecomplexHEDL");
	addBIOp<std::complex<double>,v_fesL,v_fesL,HE,tBIOp,HS_OP>("assemblecomplexHEHS");
	addBIOp<std::complex<double>,v_fesL,v_fesL,HE,tBIOpwithmass,DL_OP>("assemblecomplexHEDLwmass");
	addBIOp<std::complex<double>,v_fesL,v_fesL,HE,tBIOpwithmass,TDL_OP>("assemblecomplexHETDLwmass");

	addBIOp<std::complex<double>,v_fesL,v_fesL,HE,tBIOpcombined,SL_OP,DL_OP>("assemblecomplexHEcombinedSLDL");
	addBIOp<std::complex<double>,v_fesL,v_fesL,HE,tBIOpcombined,SL_OP,TDL_OP>("assemblecomplexHEcombinedSLTDL");
	addBIOp<std::complex<double>,v_fesL,v_fesL,HE,tBIOpcombined,HS_OP,DL_OP>("assemblecomplexHEcombinedHSDL");
	addBIOp<std::complex<double>,v_fesL,v_fesL,HE,tBIOpcombined,HS_OP,TDL_OP>("assemblecomplexHEcombinedHSTDL");

	// 2D Potential
	addPotential<std::complex<double>,v_fesL,v_fes,HE,SL_POT>("assemblecomplexHESLPot");
	addPotential<std::complex<double>,v_fesL,v_fes,HE,DL_POT>("assemblecomplexHEDLPot");

	zzzfff->Add("HMatrix", atype<HMatrixVirt<std::complex<double> > **>());
	//map_type_of_map[make_pair(atype<HMatrix<partialACA ,double>**>(), atype<double*>())] = atype<HMatrix<partialACA ,double>**>();
	map_type_of_map[make_pair(atype<HMatrixVirt<std::complex<double> >**>(), atype<Complex*>())] = atype<HMatrixVirt<std::complex<double> >**>();
}

LOADFUNC(Init_Bem)
