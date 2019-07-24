//ff-c++-LIBRARY-dep: cxx11 [mkl|blas] mpi pthread htool bemtool boost
//ff-c++-cpp-dep:
// for def  M_PI under windows in <cmath>
#define _USE_MATH_DEFINES
#include <ff++.hpp>
#include <AFunction_ext.hpp>

#include <htool/lrmat/partialACA.hpp>
#include <htool/types/matrix.hpp>
#include <htool/types/hmatrix.hpp>

// include the bemtool library .... path define in where library
//#include <bemtool/operator/block_op.hpp>  
#include <bemtool/tools.hpp>
#include <bemtool/fem/dof.hpp>
#include <bemtool/operator/operator.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>

#include "PlotStream.hpp"

extern FILE *ThePlotStream;

//#include "solve.hpp"

using namespace std;
using namespace htool;
using namespace bemtool;



template<class K>
class MyMatrix: public IMatrix<K>{
	const MeshS & ThU; // line
	const MeshS & ThV; // colunm

public:
	MyMatrix(const FESpaceS * Uh , const FESpaceS * Vh ):IMatrix<K>(Uh->Th.nv,Vh->Th.nv),ThU(Uh->Th), ThV(Vh->Th) {}

	K get_coef(const int& i, const int& j)const {return 1./(0.01+(ThU.vertices[i]-ThV.vertices[j]).norme2());}

};
/*
template<template<class> class LR, class K>
class init_Op : public E_F0mps {
public:
	Expression A;
	static const int n_name_param = 1;
	static basicAC_F0::name_and_type name_param[];
	Expression nargs[n_name_param];
	init_Op(const basicAC_F0& args, Expression param1) : A(param1) {
		args.SetNameParam(n_name_param, name_param, nargs);
	}

	AnyType operator()(Stack stack) const;
};

template<template<class> class LR, class K>
basicAC_F0::name_and_type init_Op<LR, K>::name_param[] = {
};

template<template<class> class LR, class K>
class init : public OneOperator {
public:
	init() : OneOperator(atype<HMatrix<LR,K>**>(), atype<HMatrix<LR,K>**>(), atype<Matrice_Creuse<K>*>()) { }

	E_F0* code(const basicAC_F0& args) const {
		return new init_Op<LR, K>(args, t[0]->CastTo(args[0]));
	}
};
*/



template<template<class> class LR, class K>
class assembleHMatrix : public OneOperator { public:

	class Op : public E_F0info {
	public:
		Expression a,b;

		static const int n_name_param = 8;
		static basicAC_F0::name_and_type name_param[] ;
		Expression nargs[n_name_param];
		bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
		long argl(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
		double arg(int i,Stack stack,double a) const{ return nargs[i] ? GetAny<double>( (*nargs[i])(stack) ): a;}
		KN_<long> arg(int i,Stack stack,KN_<long> a ) const{ return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
		pcommworld arg(int i,Stack stack,pcommworld a ) const{ return nargs[i] ? GetAny<pcommworld>( (*nargs[i])(stack) ): a;}
	public:
		Op(const basicAC_F0 &  args,Expression aa,Expression bb) : a(aa),b(bb) {
			args.SetNameParam(n_name_param,name_param,nargs); }
		};

		assembleHMatrix() : OneOperator(atype<const typename assembleHMatrix<LR,K>::Op *>(),
		atype<pfesS *>(),
		atype<pfesS *>()) {}

		E_F0 * code(const basicAC_F0 & args) const
		{
			return  new Op(args,t[0]->CastTo(args[0]),
			t[1]->CastTo(args[1]));
		}
	};

	template<template<class> class LR, class K>
	basicAC_F0::name_and_type  assembleHMatrix<LR, K>::Op::name_param[]= {
		{  "epsilon", &typeid(double)},
		{  "eta", &typeid(double)},
		{  "minclustersize", &typeid(long)},
		{  "maxblocksize", &typeid(long)},
		{  "mintargetdepth", &typeid(long)},
		{  "minsourcedepth", &typeid(long)},
		{		"comm", &typeid(pcommworld)},
        {        "type", &typeid(long)}
	};

	/*
	template <template<class> class LR, class K, EquationEnum Eq, BIOpKernelEnum Op, BIOpKernelEnum local_Op, int dim, typename Discretization>
	HMatrix<LR,K>* generateBIO() {
	double kappa = 1;
	Geometry node("Th.msh");
	bemtool::Mesh<dim> mesh; mesh.Load(node,0);
	Orienting(mesh);
	int nb_elt = NbElt(mesh);
	std::vector<bemtool::R3> normals= NormalTo(mesh);

	// Dof
	Dof<Discretization> dof(mesh);
	int nb_dof = NbDof(dof);
	std::vector<htool::R3> x(nb_dof);
	for (int i=0;i<nb_dof;i++){
	x[i][0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
	x[i][1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
	x[i][2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
}

BIO_Generator<BIOpKernel<Eq,Op,dim+1,Discretization,Discretization>,Discretization> generator_BIO(dof,kappa);
HMatrix<LR,K>* H = new HMatrix<LR,K>(generator_BIO,x);


//	if (apply_boundary(Op)){
//			std::vector<int> boundary=is_boundary_nodes(dof);
//			if (*max_element(boundary.begin(),boundary.end())!=0){
//					H->apply_dirichlet(boundary);
//			}
//	}

return H;
}
*/


template<template<class> class LR, class K>
AnyType SetHMatrix(Stack stack,Expression emat,Expression einter,int init)
{
	using namespace Fem2D;

	HMatrix<LR,K>** Hmat =GetAny<HMatrix<LR,K>** >((*emat)(stack));
	const typename assembleHMatrix<LR, K>::Op * mi(dynamic_cast<const typename assembleHMatrix<LR, K>::Op *>(einter));

	double epsilon=mi->arg(0,stack,htool::Parametres::epsilon);
	double eta=mi->arg(1,stack,htool::Parametres::eta);
	int minclustersize=mi->argl(2,stack,htool::Parametres::minclustersize);
	int maxblocksize=mi->argl(3,stack,htool::Parametres::eta);
	int mintargetdepth=mi->argl(4,stack,htool::Parametres::mintargetdepth);
	int minsourcedepth=mi->argl(5,stack,htool::Parametres::minsourcedepth);
	pcommworld pcomm=mi->arg(6,stack,nullptr);
    int type=mi->argl(7,stack,0);

	SetMaxBlockSize(maxblocksize);
	SetMinClusterSize(minclustersize);
	SetEpsilon(epsilon);
	SetEta(eta);
	SetMinTargetDepth(mintargetdepth);
	SetMinSourceDepth(minsourcedepth);

	




	MPI_Comm comm = pcomm ? *(MPI_Comm*)pcomm : MPI_COMM_WORLD;

	ffassert(einter);
	pfesS * pUh = GetAny< pfesS * >((* mi->a)(stack));
	FESpaceS * Uh = **pUh;
	int NUh =Uh->N;

	pfesS * pVh = GetAny< pfesS * >((* mi->b)(stack));
	FESpaceS * Vh = **pVh;
	int NVh =Vh->N;

	ffassert(Vh);
	ffassert(Uh);

	int n=Uh->NbOfDF;
	int m=Vh->NbOfDF;

	const  MeshS & ThU =Uh->Th; // line
	const  MeshS & ThV =Vh->Th; // colunm
	bool samemesh =  &Uh->Th == &Vh->Th;  // same Fem2D::Mesh

	//cout << samemesh << " " << NUh << " " << n << " " << m << endl;



	if (init)
	delete *Hmat;
    double kappa=10;
	
    //read the .msh
	std::string meshname="Th.msh";   ////// be careful
   	Geometry node(meshname);
    Mesh2D mesh;
    mesh.Load(node,0);
	Orienting(mesh);
	mesh = unbounded;
	
    
    
    
    
    
    std::cout << "Creating dof" << std::endl;
   	Dof<P1_2D> dof(mesh);
	// now the list of dof is known -> can acces to global num triangle and the local num vertice assiciated 

	
	vector<htool::R3> p1(n);
	vector<htool::R3> p2(m);
	
	
	for (int i=0; i<n; i++)
		p1[i] = {ThU.vertices[i].x, ThU.vertices[i].y, ThU.vertices[i].z};
	
    if(!samemesh){
        for (int i=0; i<m; i++)
            p2[i] = {ThV.vertices[i].x, ThV.vertices[i].y, ThV.vertices[i].z};
    }
    else
    p2=p1;
	
	// BIO_Generator
    if (type==0) {
        // equivalent myMatrix
        // MyMatrix<K> A(Uh,Vh);
        BIO_Generator<HE_SL_3D_P1xP1,P1_2D> generator_V(dof,kappa);
        *Hmat = new HMatrix<LR,K>(generator_V,p1,-1,comm);
    }
    else if (type==1) {
        
        if (mpirank==0)
        std::cout << "Loading output mesh" << std::endl;
        std::string meshname_output="Th_output.msh";   ////// be careful
        Geometry node_output(meshname_output);
        int nb_dof_output = NbNode(node_output);
      
        Potential<PotKernel<HE,SL_POT,3,P1_2D>> POT_SL(mesh,kappa);
        POT_Generator<PotKernel<HE,SL_POT,3,P1_2D>,P1_2D> generator_SL(POT_SL,dof,node_output);
        
        *Hmat = new HMatrix<LR,K>(generator_SL,p2,p1,-1,comm);
        
        
    }
	//*Hmat = generateBIO<LR,K,EquationEnum::YU,BIOpKernelEnum::HS_OP,BIOpKernelEnum::SL_OP,2,P1_2D>();

	return Hmat;
}










template<template<class> class LR, class K, int init>
AnyType SetHMatrix(Stack stack,Expression emat,Expression einter)
{ return SetHMatrix<LR,K>(stack,emat,einter,init);}

template<class R, class A, class B> R Build(A a, B b) {
	return R(a, b);
}

template<template<class> class LR, class K>
AnyType ToDense(Stack stack,Expression emat,Expression einter,int init)
{
	ffassert(einter);
	HMatrix<LR,K>** Hmat =GetAny<HMatrix<LR,K>** >((*einter)(stack));
	ffassert(Hmat && *Hmat);
	HMatrix<LR,K>& H = **Hmat;
	Matrix<K> mdense = H.to_dense_perm();
	const std::vector<K>& vdense = mdense.get_mat();

	KNM<K>* M =GetAny<KNM<K>*>((*emat)(stack));

	for (int i=0; i< mdense.nb_rows(); i++)
		for (int j=0; j< mdense.nb_cols(); j++)
			(*M)(i,j) = mdense(i,j);

	return M;
}

template<template<class> class LR, class K, int init>
AnyType ToDense(Stack stack,Expression emat,Expression einter)
{ return ToDense<LR,K>(stack,emat,einter,init);}

template<class V, template<class> class LR, class K>
class Prod {
public:
	const HMatrix<LR ,K>* h;
	const V u;
	Prod(HMatrix<LR ,K>** v, V w) : h(*v), u(w) {}

	void prod(V x) const {h->mvprod_global(*(this->u), *x);};

	static V mv(V Ax, Prod<V, LR, K> A) {
		*Ax = K();
		A.prod(Ax);
		return Ax;
	}
	static V init(V Ax, Prod<V, LR, K> A) {
		Ax->init(A.u->n);
		return mv(Ax, A);
	}

};

template<template<class> class LR, class K>
std::map<std::string, std::string>* get_infos(HMatrix<LR ,K>** const& H) {
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

template<template<class> class LR, class K>
class plotHMatrix : public OneOperator {
public:

	class Op : public E_F0info {
	public:
		Expression a;

		static const int n_name_param = 1;
		static basicAC_F0::name_and_type name_param[] ;
		Expression nargs[n_name_param];
		bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
		long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}

	public:
		Op(const basicAC_F0 &  args,Expression aa) : a(aa) {
			args.SetNameParam(n_name_param,name_param,nargs);
		}

		AnyType operator()(Stack stack) const{

			bool wait = arg(0,stack,false);

			HMatrix<LR,K>** H =GetAny<HMatrix<LR,K>** >((*a)(stack));

			PlotStream theplot(ThePlotStream);

			if (mpirank == 0 && ThePlotStream) {
				theplot.SendNewPlot();
				theplot << 3L;
				theplot <= wait;
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
				const std::vector<LR<K>*>& lrmats = (*H)->get_MyFarFieldMats();
				const std::vector<SubMatrix<K>*>& dmats = (*H)->get_MyNearFieldMats();

				int nbdense = dmats.size();
				int nblr = lrmats.size();

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

				for (int i=0;i<dmats.size();i++) {
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

				for (int i=0;i<lrmats.size();i++) {
					const LR<K>& l = *(lrmats[i]);
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

	plotHMatrix() : OneOperator(atype<long>(),atype<HMatrix<LR ,K> **>()) {}

	E_F0 * code(const basicAC_F0 & args) const
	{
		return  new Op(args,t[0]->CastTo(args[0]));
	}
};

template<template<class> class LR, class K>
basicAC_F0::name_and_type  plotHMatrix<LR, K>::Op::name_param[]= {
	{  "wait", &typeid(bool)}
};

template<template<class> class LR, class K>
void add(const char* namec) {
	Dcl_Type<HMatrix<LR ,K>**>(Initialize<HMatrix<LR,K>*>, Delete<HMatrix<LR,K>*>);
	Dcl_TypeandPtr<HMatrix<LR ,K>*>(0,0,::InitializePtr<HMatrix<LR ,K>*>,::DeletePtr<HMatrix<LR ,K>*>);
	//TheOperators->Add("<-", new init<LR,K>);
	Dcl_Type<const typename assembleHMatrix<LR, K>::Op *>();
	Add<const typename assembleHMatrix<LR, K>::Op *>("<-","(", new assembleHMatrix<LR, K>);

	TheOperators->Add("=",
	new OneOperator2_<HMatrix<LR ,K>**,HMatrix<LR ,K>**,const typename assembleHMatrix<LR, K>::Op*,E_F_StackF0F0>(SetHMatrix<LR, K, 1>));
	TheOperators->Add("<-",
	new OneOperator2_<HMatrix<LR ,K>**,HMatrix<LR ,K>**,const typename assembleHMatrix<LR, K>::Op*,E_F_StackF0F0>(SetHMatrix<LR, K, 0>));
	Global.Add(namec,"(",new assembleHMatrix<LR, K>);

	//atype<HMatrix<LR ,K>**>()->Add("(","",new OneOperator2_<string*, HMatrix<LR ,K>**, string*>(get_infos<LR,K>));
	Add<HMatrix<LR ,K>**>("infos",".",new OneOperator1_<std::map<std::string, std::string>*, HMatrix<LR ,K>**>(get_infos));

	Dcl_Type<Prod<KN<K>*, LR, K>>();
	TheOperators->Add("*", new OneOperator2<Prod<KN<K>*, LR, K>, HMatrix<LR ,K>**, KN<K>*>(Build));
	TheOperators->Add("=", new OneOperator2<KN<K>*, KN<K>*, Prod<KN<K>*, LR, K>>(Prod<KN<K>*, LR, K>::mv));
	TheOperators->Add("<-", new OneOperator2<KN<K>*, KN<K>*, Prod<KN<K>*, LR, K>>(Prod<KN<K>*, LR, K>::init));

	Global.Add("display","(",new plotHMatrix<LR, K>);

	// to dense:
	TheOperators->Add("=",
	new OneOperator2_<KNM<K>*, KNM<K>*, HMatrix<LR ,K>**,E_F_StackF0F0>(ToDense<LR, K, 1>));
	TheOperators->Add("<-",
	new OneOperator2_<KNM<K>*, KNM<K>*, HMatrix<LR ,K>**,E_F_StackF0F0>(ToDense<LR, K, 0>));
}

static void Init_Schwarz() {
	Dcl_Type<std::map<std::string, std::string>*>( );
	TheOperators->Add("<<",new OneBinaryOperator<PrintPinfos<std::map<std::string, std::string>*>>);
	Add<std::map<std::string, std::string>*>("[","",new OneOperator2_<string*, std::map<std::string, std::string>*, string*>(get_info));

	//add<partialACA,double>("assemble");
	add<partialACA,std::complex<double>>("assemblecomplex");

	zzzfff->Add("HMatrix", atype<HMatrix<partialACA ,std::complex<double> > **>());
	//map_type_of_map[make_pair(atype<HMatrix<partialACA ,double>**>(), atype<double*>())] = atype<HMatrix<partialACA ,double>**>();
	map_type_of_map[make_pair(atype<HMatrix<partialACA ,std::complex<double> >**>(), atype<Complex*>())] = atype<HMatrix<partialACA ,std::complex<double> >**>();
}

LOADFUNC(Init_Schwarz)
