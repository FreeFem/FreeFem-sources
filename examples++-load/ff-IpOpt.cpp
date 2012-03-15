/*
 *  ff-NLopt.cpp
 *  
 *
 *  Created by Sylvain Auliac on 17/01/12.
 *
 */
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
//ff-c++-LIBRARY-dep:  Ipopt mumps-seq blas  libseq  fc  


#include  <iostream>
#include <stack>
#include <vector>
using namespace std;
#include "ff++.hpp"


#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include <list>
#include <set>
#include <map>
#include <cstdarg>


extern Block *currentblock;

typedef double R;
typedef KN_<R> Rn_;
typedef KN<R> Rn;
typedef KNM_<R> Rnm_;
typedef KNM<R> Rnm;





/*****************************************************************************************************************************
 *	Some misc. function usefull later...
 *****************************************************************************************************************************/
 
//A variadic function to add an undefinite number of elements to a set of short int
//This is used to define the set of named parameter which are not used when certain assumptions
//upon the optimization poblem functions are met 
void AddElements(std::set<unsigned short> &_set,int amount,int first,...)
{ 
	int elem=0;
  va_list vl;
  va_start(vl,first);
	_set.insert(first);
  for (int i=1;i<amount;i++)
  {
    elem=va_arg(vl,int);
		_set.insert(elem);
		
  }
  va_end(vl);

}
//A raw pointer cleaner
template<class T> inline void clean(T *p) {if(p) {delete p; p=0;} }
inline bool operator<=(const std::pair<int,int> &l,const std::pair<int,int> &r)
{
	return (l.first < r.first) || (l.first==r.first && l.second <= r.second);
}
inline bool XOR(bool a,bool b) {return (!a && b) || (a && !b);}
inline bool NXOR(bool a,bool b) {return !XOR(a,b);}
inline void Sonde(int i) {cout << "sonde " << i << endl;}



template<class K> class ffcalfunc  //   to call the freefem function .. J, constraints, and derivatives
{ 
	public:
		Stack stack;
		Expression JJ,theparame;
		ffcalfunc(const ffcalfunc &f) : stack(f.stack),JJ(f.JJ),theparame(f.theparame) {}
		ffcalfunc(Stack s,Expression JJJ,Expression epar) : stack(s),JJ(JJJ), theparame(epar) {}
		K J(Rn_  x) const 
		{
			KN<double> *p=GetAny<KN<double> *>( (*theparame)(stack) );
			*p=x;
			K ret= GetAny<K>( (*JJ)(stack));
			//cout << "call to ffcalfunc.J with " << *p << " and ret=" << ret << endl;
			WhereStackOfPtr2Free(stack)->clean();
			return  ret; 
		}
};
template<> class ffcalfunc<Matrice_Creuse<R> *>
{
	public:
		typedef Matrice_Creuse<R> *K;
		Stack stack;
		ffcalfunc(const ffcalfunc &f) : stack(f.stack) {}
		ffcalfunc(Stack s) : stack(s) {}
		virtual K J(Rn_) const = 0;
		virtual K J(Rn_,double,Rn_) const = 0;
		virtual bool NLCHPEnabled() const = 0; //Non Linear Constraints Hessian Prototype 
};


class GeneralSparseMatFunc : public ffcalfunc<Matrice_Creuse<R> *>
{
	private:
		typedef ffcalfunc<Matrice_Creuse<R> *> FFF;
	public:
		Expression JJ,param,paramlm,paramof;
		GeneralSparseMatFunc(const GeneralSparseMatFunc &f) : FFF(f),JJ(f.JJ),param(f.param),paramlm(f.paramlm),paramof(f.paramof) {};
		GeneralSparseMatFunc(Stack s,Expression JJJ,Expression epar,Expression eparof=0,Expression eparlm=0) 
			: FFF(s),JJ(JJJ),param(epar),paramlm(eparlm),paramof(eparof)
		{ffassert(NXOR(paramlm,paramof));}
		bool NLCHPEnabled() const {return paramlm && paramof;}
		K J(Rn_  x) const 
		{
			KN<double> *p=GetAny<KN<double> *>( (*param)(stack) );
			*p=x;
			K ret= GetAny<K>( (*JJ)(stack));
			//cout << "call to ffcalfunc.J with " << *p << " and ret=" << ret << endl;
			WhereStackOfPtr2Free(stack)->clean();
			return  ret; 
		}
		K J(Rn_  x,double of,Rn_ lm) const 
		{
			if(paramlm && paramof)
			{
				KN<double> *p=GetAny<KN<double> *>( (*param)(stack) );
				double *pof=GetAny<double *>( (*paramof)(stack) );
				KN<double > *plm=GetAny<KN<double> *>( (*paramlm)(stack) );
				*p=x;
				*pof=of;
				*plm=lm;
				K ret= GetAny<K>( (*JJ)(stack));
				//cout << "call to ffcalfunc.J with " << *p << " and ret=" << ret << endl;
				WhereStackOfPtr2Free(stack)->clean();
				return  ret; 
			}
			else return J(x);
		}
};

class ConstantSparseMatFunc : public ffcalfunc<Matrice_Creuse<R> *>
{
	private:
		typedef ffcalfunc<Matrice_Creuse<R> *> FFF;
	public:
		Expression M; //Expression of the matrix
		ConstantSparseMatFunc(const ConstantSparseMatFunc &f) : FFF(f),M(f.M) {}
		ConstantSparseMatFunc(Stack s,Expression _M) : FFF(s),M(_M) {}
		bool NLCHPEnabled() const {return false;}
		K J(Rn_) const
		{
			K ret = GetAny<K>( (*M)(stack) );
			WhereStackOfPtr2Free(stack)->clean();
			return ret;
		}
		K J(Rn_ x,double,Rn_) const {return J(x);}
};
		





typedef ffcalfunc<double> ScalarFunc;
typedef ffcalfunc<Rn> VectorFunc;
typedef ffcalfunc<Rnm> FullMatrixFunc;
typedef ffcalfunc<Matrice_Creuse<R>* > SparseMatFunc;


class SparseMatStructure
{
	public:
		typedef std::pair<int,int> Z2;
		typedef std::set<Z2> Structure;
		typedef std::pair<KN<int>,KN<int> > Zn2;
		typedef Structure::const_iterator const_iterator;
		typedef Structure::iterator iterator;
		
		SparseMatStructure(bool _sym=0) : structure(),sym(_sym),n(0),m(0),raws(0),cols(0) {}
		SparseMatStructure(Matrice_Creuse<R> const * const M,bool _sym=0) : structure(),sym(_sym),n(M->N()),m(M->M()),raws(0),cols(0) {this->AddMatrix(M);}
		template<class INT> SparseMatStructure(const KN<INT> &I,const KN<INT> &J,bool _sym=0) : structure(),sym(_sym),n(I.max()),m(J.max()),raws(0),cols(0) {this->AddArrays(I,J);}
		~SparseMatStructure() {if(raws) delete raws; if(cols) delete cols;}
		
		const_iterator begin() const {return structure.begin();}
		iterator begin() {return structure.begin();}
		const_iterator end() const {return structure.end();}
		iterator end() {return structure.end();}
		//Structure& operator()() {return structure;}
		//const Structure& operator()() const {return structure;}
		bool empty() const {return structure.empty() && !raws && !cols;}
		int N() const {return n;}
		int M() const {return m;}
		
		SparseMatStructure& clear() {structure.clear(); if(raws) delete raws; if(cols) delete cols; sym=false; n=0; m=0; return *this;}
		int size() const {return structure.size() ? structure.size() : (raws ? raws->N() : 0);}
		SparseMatStructure& AddMatrix(Matrice_Creuse<R> const * const);
		template<class INT> SparseMatStructure& AddArrays(const KN<INT> &,const KN<INT> &);
		SparseMatStructure& ToKn(bool emptystruct=false);
		
		
		KN<int> & Raws() {return *raws;}
		KN<int> const & Raws() const {return *raws;}
		KN<int> & Cols() {return *cols;}
		KN<int> const & Cols() const {return *cols;}
		
	private:
		int n,m;
		Structure structure;
		bool sym;
		//Zn2 *array_structure;
		KN<int> *raws,*cols;
};

SparseMatStructure& SparseMatStructure::ToKn(bool emptystruct)
{
	if(raws) delete raws;
	if(cols) delete cols;
	raws = new KN<int>(structure.size());
	cols = new KN<int>(structure.size());
	int k=0;
	for(const_iterator i=begin();i!=end();++i) {(*raws)[k]=i->first; (*cols)[k]=i->second; ++k;}
	if(emptystruct) structure.clear();
	return *this;
}

SparseMatStructure& SparseMatStructure::AddMatrix(Matrice_Creuse<R> const * const _M)
{
	n = n > _M->N() ? n : _M->N();
	m = m > _M->M() ? m : _M->M();
	MatriceMorse<R> const * const M = dynamic_cast<MatriceMorse<R> const * const> (&(*_M->A));
	if(!sym || (sym && M->symetrique))
	{
		for(int i=0;i < M->N;++i)
		{
			for(int k=M->lg[i]; k < M->lg[i+1]; ++k) structure.insert(Z2(i,M->cl[k]));
		}
	}
	else // sym && !M->symetrique
	{
		for(int i=0;i<M->N;++i)
		{
			for(int k=M->lg[i]; k < M->lg[i+1]; ++k) if(i >= M->cl[k]) structure.insert(Z2(i,M->cl[k]));
		}
	}
	return *this;
}
template<class INT> SparseMatStructure& SparseMatStructure::AddArrays(const KN<INT> &I,const KN<INT> &J)
{
	ffassert(I.N()==J.N());
	n = n > I.max()+1 ? n : I.max()+1;
	m = m > J.max()+1 ? m : J.max()+1;
	if(!sym) for(int k=0;k<I.N();++k) structure.insert(Z2(I[k],J[k]));
	else for(int k=0;k<I.N();++k) if(I[k]>=J[k]) structure.insert(Z2(I[k],J[k]));
	return *this;
}

using namespace Ipopt;

class ffNLP : public TNLP
{
	public:
		ffNLP() : xstart(0) {}
		ffNLP(Rn &,const Rn &,const Rn &,const Rn &,const Rn &,ScalarFunc*, VectorFunc*, SparseMatFunc*, VectorFunc*, SparseMatFunc*);
		ffNLP(Rn &,const Rn &,const Rn &,const Rn &,const Rn &,ScalarFunc*, VectorFunc*, SparseMatFunc*, VectorFunc*, SparseMatFunc*, int ,int ,int);
		virtual ~ffNLP();
		
		bool get_nlp_info(Index&, Index&, Index&, Index&, IndexStyleEnum&); //the IPOPT methods
		bool get_bounds_info(Index, Number*, Number*, Index, Number*, Number*);
		bool get_starting_point(Index, bool, Number*,bool , Number* , Number*,Index , bool ,Number* );
		bool eval_f(Index, const Number*, bool, Number&);
		bool eval_grad_f(Index, const Number*, bool, Number*);
		bool eval_g(Index, const Number*, bool, Index, Number*);
		bool eval_jac_g(Index, const Number*, bool,Index, Index, Index*, Index *,Number*);
		bool eval_h(Index, const Number*, bool ,Number , Index , const Number*,bool, Index, Index*,Index*, Number*);
		void finalize_solution(SolverReturn, Index, const Number*, const Number*, const Number*, Index, const Number*, const Number*, Number,
				 const IpoptData* ip_data,
				 IpoptCalculatedQuantities* ip_cq);

		template<class INT> ffNLP& SetHessianStructure(const KN<INT> &,const KN<INT> &,bool reset=0);
		template<class INT> ffNLP& SetJacobianStructure(const KN<INT> &,const KN<INT> &,bool reset=0);
		enum Level {do_nothing,user_defined, one_evaluation, basis_analysis};
		ffNLP& BuildMatrixStructures(Level,Level,int);
		ffNLP& EnableCheckStruct() {checkstruct=true; return *this;}
		ffNLP& DisableCheckStruct() {checkstruct=false; return *this;}
		
		Rn lambda_start,x_start;
		double sigma_start;
		
	private:
		//algorithm datas
		Rn *xstart,xl,xu,gl,gu;
		ScalarFunc *fitness;
		VectorFunc *dfitness,*constraints;
		SparseMatFunc *hessian,*dconstraints;
		int mm,nnz_jac,nnz_h;
		double final_value;
		//bool sym;
		bool checkstruct;
		SparseMatStructure HesStruct,JacStruct;
		
		//some static functions...
		template<class A,class B> static void KnToPtr(const KN<A> &a,B *b) {for(int i=0;i<a.N();++i) b[i]=a[i];}
		template<class A,class B> static void KnFromPtr(KN<A> &a,B const *b) {for(int i=0;i<a.N();++i) a[i]=b[i];}
		static int FindIndex(const KN<int> &irow,const KN<int> & jrow,int i,int j,int kmin,int kmax);
};


ffNLP::ffNLP(Rn &x,const Rn &_xl,const Rn &_xu,const Rn &_gl,const Rn &_gu,ScalarFunc * _fitness,VectorFunc * _dfitness,SparseMatFunc * _hessian, 
						 VectorFunc * _constraints,SparseMatFunc * _dconstraints) : 
				  xstart(&x), xl(_xl), xu(_xu), gl(_gl), gu(_gu),final_value(299792458.),//sym(0),unsymind(),
					fitness(_fitness), dfitness(_dfitness), constraints(_constraints),
					hessian(_hessian), dconstraints(_dconstraints),mm(-1),nnz_jac(-1),nnz_h(-1),
					HesStruct(true),JacStruct(false),sigma_start(1.),lambda_start(),x_start(x),checkstruct(1) {}


ffNLP::ffNLP(Rn &x,const Rn &_xl,const Rn &_xu,const Rn &_gl,const Rn &_gu,ScalarFunc * _fitness,VectorFunc * _dfitness,SparseMatFunc * _hessian,
						 VectorFunc * _constraints,SparseMatFunc * _dconstraints, int _mm,int _nnz_jac,int _nnz_h) : 
				  xstart(&x), xl(_xl), xu(_xu), gl(_gl), gu(_gu),hessian(_hessian),final_value(299792458.),//sym(0),unsymind(),
					fitness(_fitness),dfitness(dfitness),constraints(_constraints),dconstraints(_dconstraints),
					mm(_mm),nnz_jac(_nnz_jac),nnz_h(_nnz_h),HesStruct(true),JacStruct(false),sigma_start(1.),lambda_start(),x_start(x),checkstruct(1) {}

ffNLP::~ffNLP()
{
	/*
	clean(fitness);
	clean(dfitness);
	clean(constraints);
	clean(hessian);
	clean(dconstraints);
	*/
}

template<class INT> ffNLP& ffNLP::SetHessianStructure(const KN<INT> &I,const KN<INT> &J,bool reset)
{
	if(reset) HesStruct.clear();
	HesStruct.AddArrays(I,J);
	return *this;
}
template<class INT> ffNLP& ffNLP::SetJacobianStructure(const KN<INT> &I,const KN<INT> &J,bool reset)
{
	if(reset) JacStruct.clear();
	JacStruct.AddArrays(I,J);
	return *this;
}
ffNLP& ffNLP::BuildMatrixStructures(Level hlvl, Level jlvl,int _mm)
{
	if(jlvl!=do_nothing && dconstraints)
	{
		if(jlvl==user_defined) ffassert(JacStruct.size());
		else if((jlvl==one_evaluation || jlvl==basis_analysis) && dconstraints) JacStruct.AddMatrix(dconstraints->J(x_start));
	}
	if(hlvl!=do_nothing && hessian)
	{
		if(hlvl==user_defined) ffassert(HesStruct.size());
		else if(hlvl==one_evaluation || !hessian->NLCHPEnabled() ) HesStruct.AddMatrix(hessian->J(x_start,sigma_start,lambda_start));
		else if(hlvl==basis_analysis)
		{
			{
				Rn lambda(_mm,0.);
				HesStruct.AddMatrix(hessian->J(x_start,1.,lambda));
			}
			for(int i=0;i<_mm;++i)
			{
				Rn lambda(_mm,0.);
				lambda[i] = 1.;
				HesStruct.AddMatrix(hessian->J(x_start,0.,lambda));
				lambda[i] = 0.;
			}
		}
	}
	JacStruct.ToKn();
	HesStruct.ToKn();
	return *this;
}
int ffNLP::FindIndex(const KN<int> &irow,const KN<int> &jcol,int i,int j,int kmin,int kmax)
{
	//cout << "Trying to find (" << i << ',' << j << ") in :" << irow << jcol << " - kmin=" << kmin << " and kmax=" << kmax << endl;
	typedef std::pair<int,int> Z2;
	Z2 ij(i,j),ijmin(irow[kmin],jcol[kmin]),ijmax(irow[kmax],jcol[kmax]);
	if(abs(kmin-kmax)<=1)
	{
		if(ij==ijmin) return kmin;
		else if(ij==ijmax) return kmax;
		else return -1;
	}
	else
	{
		int knew = (kmin + kmax) / 2;
		Z2 ijnew(irow[knew],jcol[knew]);
		if(ij <= ijnew) return FindIndex(irow,jcol,i,j,kmin,knew);
		else return FindIndex(irow,jcol,i,j,knew,kmax);
	}
}


bool ffNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	bool ret=true;
	n = xstart ? xstart->N() : (ret=0);
	//set(m,mm,constraints,xstart,ret);
	//set(nnz_jac_g,nnz_jac,dconstraints,xstart,ret);
	//set(nnz_h_lag,nnz_h,hessian,xstart,ret);
	//if(JacStruct.empty() && constraints) BuildMatrixStructures(do_nothing,one_evaluation);
	//if(HesStruct.empty()) BuildMatrixStructures(one_evaluation,do_nothing);
	mm = m = constraints ? JacStruct.N() : 0;
	nnz_jac = nnz_jac_g = constraints ? JacStruct.size() : 0;
	nnz_h = nnz_h_lag = HesStruct.size();
	index_style = TNLP::C_STYLE;
	return ret;
}

bool ffNLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
	//cout << "n=" << n << " m=" << m << " mm=" << mm << " g_l.N()=" << gl.N() << " g_u.N()=" << gu.N() << endl;
	//assert(gl.N()==mm);
	//assert(gu.N()==mm);
	KnToPtr(xl,x_l);
	KnToPtr(xu,x_u);
	if(mm) KnToPtr(gl,g_l);
	if(mm) KnToPtr(gu,g_u);
	/* DEBUG
	cout << "constraints lower bound = (";
	for(int i=0;i<m;++i) cout << g_l[i] <<  (i<m-1 ? ',':')');
	cout << endl << "constraints upper bound = (";
	for(int i=0;i<m;++i) cout << g_u[i] <<  (i<m-1 ? ',':')');
	cout << endl;*/
	return true;
}
bool ffNLP::get_starting_point(Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U,Index m, bool init_lambda,Number* lambda)
{
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);
	assert(xstart->N() == n);
	KnToPtr(*xstart,x);
	return true;
}
bool ffNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
	assert(n == xstart->N());
	Rn X(n);
	KnFromPtr(X,x);
	obj_value = fitness->J(X);
	return true;
}
bool ffNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == xstart->N());
	Rn X(n);
	KnFromPtr(X,x);
	Rn _grad_f=dfitness->J(X);
	KnToPtr(_grad_f,grad_f);
  return true;
}
bool ffNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
	Rn X(n);
	KnFromPtr(X,x);
	if(constraints)
	{
		Rn _g=constraints->J(X);
		KnToPtr(_g,g);
	}
	return true;
}
bool ffNLP::eval_jac_g(Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values)
{
	assert(n==xstart->N());
	Rn X(n);
	if(x) KnFromPtr(X,x); else X=*xstart;
	
	if(values==0)
	{
		int k=0;
		for(SparseMatStructure::const_iterator i=JacStruct.begin(); i != JacStruct.end(); ++i)
		{
			iRow[k] = i->first;
			jCol[k] = i->second;
			++k;
		}
	}
	else 
	{
		Matrice_Creuse<R>* M = dconstraints->J(X);
		MatriceMorse<R> *MM = dynamic_cast<MatriceMorse<R>* >(&(*M->A));  //ugly!
		for(int i=0;i<MM->N;++i)
		{
			for(int k=MM->lg[i]; k < MM->lg[i+1]; ++k)
			{
				if(checkstruct)
				{
					int kipopt = FindIndex(JacStruct.Raws(),JacStruct.Cols(),i,MM->cl[k],0,nele_jac-1);
					if(kipopt>=0) values[kipopt] = MM->a[k];
				}
				else values[k] = MM->a[k];
			}
		}
	}
	return true;
}


bool ffNLP::eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda, Index nele_hess, Index* iRow,Index* jCol, Number* values)
{
	Rn X(n),L(m);
	if(x) KnFromPtr(X,x); else X=*xstart;
	if(lambda) KnFromPtr(L,lambda); else L=0.;
	
	bool NLCHPE = hessian->NLCHPEnabled();
	Number _obj_factor = NLCHPE ? 1. : obj_factor;
	if(values==0)
	{
		int k=0;
		for(SparseMatStructure::const_iterator i=HesStruct.begin(); i != HesStruct.end(); ++i)
		{
			iRow[k] = i->first;
			jCol[k] = i->second;
			++k;
		}
	}
	else
	{
		Matrice_Creuse<R>* M=0;
		if(NLCHPE) M=hessian->J(X,obj_factor,L); else M=hessian->J(X);
		MatriceMorse<R> *MM = dynamic_cast<MatriceMorse<R>* >(&(*M->A));//ugly!
		if(checkstruct)
		{
			for(int i=0;i<MM->N;++i)
			{
				for(int k=MM->lg[i]; k < MM->lg[i+1]; ++k)
				{
					int kipopt = FindIndex(HesStruct.Raws(),HesStruct.Cols(),i,MM->cl[k],0,nele_hess-1);
					if(kipopt>=0) values[kipopt] = _obj_factor * (MM->a[k]);
					//else values[k] = (hessian->paramof &&hessian->paramlm ? 1. : obj_factor) * (MM->a[k]);
				}
			}
		}
		else if(! MM->symetrique)
		{
			for(int i=0,kipopt=0;i<MM->N;++i)
			{
				for(int k=MM->lg[i]; k < MM->lg[i+1]; ++k)
				{
					if(i >= MM->cl[k]) 
					{
						values[kipopt] = _obj_factor * (MM->a[k]);
						++kipopt;
					}
				}
			}
		}
		else
		{
			for(int i=0;i<MM->N;++i)
			{
				for(int k=MM->lg[i]; k < MM->lg[i+1]; ++k) values[k] = _obj_factor * (MM->a[k]);
			}
		}
	}

  return true;
}


void ffNLP::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,const IpoptData* ip_data,
				 IpoptCalculatedQuantities* ip_cq)
{
	KnFromPtr(*xstart,x);
	KnFromPtr(lambda_start,lambda);
	
	final_value = obj_value;
}


static ffNLP::Level ToLevel(long i)
{
	switch(i){
		case 0:
			return ffNLP::user_defined;
			break;
		case 1:
			return ffNLP::one_evaluation;
			break;
		case 2:
			return ffNLP::basis_analysis;
			break;
		default:
			return ffNLP::do_nothing;
			break;
		}
}
inline void SONDE() {static int i=1; cout << "SONDE " << i << endl; ++i;}


//General case : (no_assumption)
//  -handles all possibilities
//  -fitness function, its gradient, the lagrangian hessian, constraint and its jacobian are all passed with functions
//  -accepts two different prototype for the lagrangian hessian (real[int]&,real,real[int]&) for non linear constraints
//   as well as (real[int] &) which is meant to be used when the constraints are affine (no contribution to the hessian
//   since second order derivatives of constraints are null)

//Affine constraints case: (affine_g)
//  -handles the case when constraints have null second order derivatives
//  -fitness function, its gradient and hessian, and constraints are still passed by functions
//  -as it is constant, constraints jacobian can be directly passed with a matrix

//Constant lagrangian hessian case : (quadratic_f_affine_g)
//  -handles the case of a constant lagrangian hessian (which necessary means that constraints are affine, and fitness function is quadratic)
//  -only fitness function, gradient and constraints are passed by functions
//  -all matricial functions are directly passed with a constant matrix

//Constant lagrangian with non affine constraints  : doesn't exist --> spurious case

enum Assumption {no_assumption, affine_g, quadratic_f_affine_g, spurious_assumption, bfgs, bfgs_affine_g};
const bool with_constraints=true,without_constraints=false;
template<Assumption A,bool WC> struct Case
{
	Case() {}
	static const Assumption a=A;
	static const bool wc=WC;
};




class OptimIpopt : public OneOperator
{
	public:
		const Assumption A;
		const bool WC;
		
		class E_Ipopt : public E_F0mps
		{
			private:
				bool first;
				std::set<unsigned short> unused_name_param;
				void InitUNP();
				bool CompletelyNonLinearConstraints;
			public:
				const Assumption A;
				const bool WC;
				static basicAC_F0::name_and_type name_param[];
				static const int n_name_param=12;
				Expression nargs[n_name_param];
				Expression X;
				mutable Rn lm;
				C_F0 L_m;
				C_F0 inittheparam,theparam,closetheparam;
				C_F0 initobjfact,objfact;
				Expression JJ,GradJ,Constraints,GradConstraints,Hessian;
				bool arg(int i,Stack stack,bool a) const {return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
				long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
				R arg(int i,Stack stack,R a) const{ return nargs[i] ? GetAny<R>( (*nargs[i])(stack) ): a;}
				Rn_ arg(int i,Stack stack,Rn_ a) const {return nargs[i] ? GetAny<Rn_>((*nargs[i])(stack)) : a;}
				template<typename T> T Arg(int i,Stack s) const {return GetAny<T>( (*nargs[i])(s));}
				
				E_Ipopt(const basicAC_F0 & args,Assumption a,bool wc) : CompletelyNonLinearConstraints(true),lm(),L_m(CPValue(lm)),A(a),WC(wc),first(true),unused_name_param()
				{
					if(first) {InitUNP(); first=false;}
					if(A==spurious_assumption && WC) CompileError("IPOPT: impossible to have a constant lagrangian hessian and non linear constraints at the same time (see documentation)\n or if your constraints have a constant jacobian, you should directly pass the matrix to IPOPT, just as you did for the hessian.");
					int nbj= args.size()-1;
					Block::open(currentblock); // make a new block to 
					X = to<Rn*>(args[nbj]);
					C_F0 X_n(args[nbj],"n");
					//  the expression to init the theparam of all 
					inittheparam = currentblock->NewVar<LocalVariable>("the parameter",atype<KN<R> *>(),X_n);
					initobjfact = currentblock->NewVar<LocalVariable>("objective factor",atype<double *>());
					theparam = currentblock->Find("the parameter"); //  the expression for the parameter
					objfact = currentblock->Find("objective factor");
					args.SetNameParam(n_name_param,name_param,nargs);
					const  Polymorphic * opJ=0,*opdJ=0,*opH=0,*opG=0,*opjG=0;
					if (nbj>0)
					{
						opJ = dynamic_cast<const  Polymorphic *>(args[0].LeftValue());
						assert(opJ);
						opdJ= dynamic_cast<const Polymorphic *> (args[1].LeftValue());
						assert(opdJ);
						if(A!=quadratic_f_affine_g && A!=bfgs && A!=bfgs_affine_g)
						{
							opH = dynamic_cast<const Polymorphic *> (args[2].LeftValue());
							assert(opH);
						}
						if(WC)
						{
							opG = dynamic_cast<const Polymorphic *> (args[nbj-2].LeftValue());
							assert(opG);
							if(A!=affine_g && A!=quadratic_f_affine_g && A!=bfgs_affine_g)
							{
								opjG= dynamic_cast<const Polymorphic *> (args[nbj-1].LeftValue());
								assert(opjG); 
							}
						}
					}      
					ArrayOfaType hprototype2(atype<KN<R> *>(),atype<double>(),atype<KN<R>*>()),hprototype1(atype<KN<R> *>());
					JJ =	to<R>(C_F0(opJ,"(",theparam));
					GradJ = to<Rn_>(C_F0(opdJ,"(",theparam));
					if(A!=quadratic_f_affine_g && A!=bfgs && A!=bfgs_affine_g)
					{
						if(opH->Find("(",hprototype2))
						{
							CompletelyNonLinearConstraints = true;
							Hessian = to<Matrice_Creuse<R>* >(C_F0(opH,"(",theparam,objfact,L_m));
						}
						else if(opH->Find("(",hprototype1))
						{
							CompletelyNonLinearConstraints = false; //When constraints are affine, lagrange multipliers are not used in the hessian, obj_factor is also hidden to the user
							Hessian = to<Matrice_Creuse<R>* >(C_F0(opH,"(",theparam));
						}
						else CompileError("Error, wrong hessian function prototype. Must be either (real[int] &) or (real[int] &,real,real[int] &)");
					}
					else if(A!=bfgs && A!=bfgs_affine_g) Hessian = to<Matrice_Creuse<R> *>(args[2]);
					if(WC)
					{
						Constraints = to<Rn_>(C_F0(opG,"(",theparam));
						if(A==affine_g || A==quadratic_f_affine_g || A==bfgs_affine_g) GradConstraints = to<Matrice_Creuse<R> *>(args[nbj-1]);
						else GradConstraints = to<Matrice_Creuse<R>*>(C_F0(opjG,"(",theparam));
					}
					closetheparam=currentblock->close(currentblock);   // the cleanning block expression 
				}
				
				
				virtual AnyType operator()(Stack stack)  const
				{
					//cout << "ASSUMPTION = " << A << "  -  WC = " << WC << endl;
					double cost = 299792458.;
					WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH mars 2005 
					Rn &x = *GetAny<Rn *>((*X)(stack));	
					long n=x.N();	
					bool warned=false;
					for(int i=0;i<n_name_param;++i) 
						if(nargs[i] && unused_name_param.find(i)!=unused_name_param.end())
						{
							cout << "IPOPT Warning: named parameter " << name_param[i].name << " is useless for the problem you set (indications may appear further...)." << endl;
							warned = true;
						}
					if(warned)
					{
						if(!WC && (nargs[2] || nargs[3])) cout << "  ==> Some constraints bounds have been defined while no constraints function has been passed." << endl;
						if(!WC && nargs[4]) cout << "  ==> You provided a structure for the constraints jacobian but there is no constraint function." << endl;
						if(!WC && nargs[6]) cout << "  ==> Unconstrained problem make the use of " << name_param[6].name << " pointless (see the documentation for more details)" << endl;
						if(!WC && A==bfgs && nargs[8])
						{
							cout << "  ==> " << name_param[8].name << " is useless because there should not be any function returning matrix in your problem," << endl;
							cout << "      (2 functions can only be J and dJ). You may as well have forgotten one function (ipopt will certainly crash if so)." << endl;
						}
						if(WC && A==bfgs_affine_g && nargs[8])
						{
							cout << "  ==> " << name_param[8].name << " is useless because there should not provide any non constant matrix in your problem," << endl;
							cout << "      (3 functions can only be J, dJ and constraints). You may as well have forgotten one function (ipopt will certainly crash if so)." << endl;
						}
						if((A==affine_g || A==quadratic_f_affine_g || A==bfgs_affine_g) && WC && nargs[4])
						{
							cout << "  ==> your constraints jacobian is a constant matrix, thus there is no need to specify its structure with " << name_param[4].name << endl;
							cout << "      since it is contained in the matrix you passed." << endl;
						}
						/*if((A==affine_g || A==quadratic_f_affine_g) && WC && nargs[6])
						{
							cout << "  ==> " << name_param[6].name << " is meant to be used for the automatic determination of the lagrangian hessian structure," << endl;
							cout << "      but all of your constraints have null hessian, thus they won't contribute to it (see the documentation for more details)." << endl;
						}*/
						if(A==quadratic_f_affine_g && nargs[5])
						{
							cout << "  ==> your lagrangian hessian is a constant matrix, thus there is no need to specify its structure with " << name_param[5].name << endl;
							cout << "      since it is contained in the matrix you passed." << endl;
						}
						if(A==quadratic_f_affine_g && nargs[7])
						{
							cout << "  ==> " << name_param[7].name << " will be ignored since all matricial objects are constants and constraints do not" << endl;
							cout << "      contribute to the hessian, matrix structure determination is trivial." << endl;
						}
						if(A==bfgs_affine_g && (nargs[7] || nargs[8]))
						{
							cout << "  ==> " << name_param[7].name << "or" << name_param[7].name << " will be ignored since the only matrix object you have passed is constant " << endl;
							cout << "      the structure determination is contained in it." << endl;
						}
						if(A==quadratic_f_affine_g && nargs[8])
						{
							//if(Arg<bool>(8,stack))
							cout << "  ==> no need to use " << name_param[8].name << " since all matricial objects are constants, structures won't change through the algorithm," << endl;
							cout << "      it is automatically set to the default disabling value." << endl;
						}
						if(A==bfgs && nargs[5])
						{
							cout << "  ==> the number of functions you passed to ipopt is such that the LBFGS mode has been enabled, thus making " << endl;
							cout << "      " << name_param[5].name << " useless. You may also have forgoten a function (ipopt will certainly crash if so)." << endl;
						}
					}
					
					
					long iprint = verbosity;	
					ScalarFunc ffJ(stack,JJ,theparam);
					VectorFunc ffdJ(stack,GradJ,theparam);
					SparseMatFunc * ffH=0;
					if(CompletelyNonLinearConstraints) ffH = new GeneralSparseMatFunc(stack,Hessian,theparam,objfact,L_m);
					else ffH = new GeneralSparseMatFunc(stack,Hessian,theparam);
					VectorFunc *ffC = WC ? new ffcalfunc<Rn>(stack,Constraints,theparam) : 0;
					SparseMatFunc *ffdC = WC ?  new GeneralSparseMatFunc(stack,GradConstraints,theparam) : 0; 
					
					Rn xl(n),xu(n),gl(nargs[2] ? Arg<Rn_>(2,stack).N() : 0),gu(nargs[3] ? Arg<Rn_>(3,stack).N() : 0);
					if(WC && (gl.N()+gu.N())==0) cout << "IPOPT Warning : constrained problem without constraints bounds" << endl;
					int mmm=gl.N()>gu.N() ? gl.N() : gu.N();
					Rn_ *lag_mul=0;//Rn(mmm,1.);
					//int niter=arg(6,stack,100L);
					int autostructmode = (A==quadratic_f_affine_g || A==bfgs_affine_g )? ffNLP::one_evaluation : arg(7,stack,1L);
					bool checkindex = ((A==quadratic_f_affine_g || A==bfgs_affine_g) ? false : arg(8,stack,false)), cberror=false;
					
					
					if(nargs[0]) xl=Arg<Rn_>(0,stack); else xl=-1.e19;
					if(nargs[1]) xu=Arg<Rn_>(1,stack); else xu=1.e19;
					if(nargs[2]) gl=Arg<Rn_>(2,stack); else {gl.resize(mmm); gl=-1.e19;}
					if(nargs[3]) gu=Arg<Rn_>(3,stack); else {gu.resize(mmm); gu=1.e19;}
					const E_Array * ejacstruct = (WC && A==no_assumption && nargs[4]) ? dynamic_cast<const E_Array *> (nargs[4]) : 0,
												* ehesstruct = (A!=bfgs && A!=bfgs_affine_g && A!=quadratic_f_affine_g && nargs[5]) ? dynamic_cast<const E_Array *> (nargs[5]) : 0;
					//Expression LM = 0;
					if(nargs[6] && WC)
					{
						lag_mul = new Rn_(GetAny<Rn_>((*nargs[6])(stack)));
					}
					
					SmartPtr<TNLP> optim = new ffNLP(x,xl,xu,gl,gu,&ffJ,&ffdJ,ffH,ffC,ffdC);
					ffNLP * _optim = dynamic_cast<ffNLP *> (&(*optim));
					assert(_optim);
					if(WC && nargs[6]) _optim->lambda_start = *lag_mul;
					else if(WC)
					{
						_optim->lambda_start.resize(mmm);
						_optim->lambda_start = 1.;
					}
					_optim->sigma_start = 1.;
					
					if(ejacstruct)
					{
						if(ejacstruct->nbitem()!=2) ExecError("\nSorry, we were expecting an array with two componants in structjac=[iraw,jcol]");
						if((*ejacstruct)[0].left() != atype<KN<long> *>()) CompileError("Sorry, array componants in structjac=[iraw,jcol] must be integral type vectors");
						if((*ejacstruct)[1].left() != atype<KN<long> *>()) CompileError("Sorry, array componants in structjac=[iraw,jcol] must be integral type vectors");
						Expression raws = (*ejacstruct)[0], cols = (*ejacstruct)[1];
						_optim->SetJacobianStructure(*GetAny<KN<long>*>((*raws)(stack)),*GetAny<KN<long>*>((*cols)(stack)),true);
					}
					if(ehesstruct)
					{
						if(ehesstruct->nbitem()!=2) ExecError("\nSorry, we were expecting an array with two componants in structhess=[iraw,jcol]");
						if((*ehesstruct)[0].left() != atype<KN<long> *>()) CompileError("Sorry, array componants in structhess=[iraw,jcol] must be integral type vectors");
						if((*ehesstruct)[1].left() != atype<KN<long> *>()) CompileError("Sorry, array componants in structhess=[iraw,jcol] must be integral type vectors");
						Expression raws = (*ehesstruct)[0], cols = (*ehesstruct)[1];
						_optim->SetHessianStructure(*GetAny<KN<long>*>((*raws)(stack)),*GetAny<KN<long>*>((*cols)(stack)),true);
					}
					ffNLP::Level lh=ehesstruct ? ffNLP::user_defined : ToLevel(autostructmode),lj=ejacstruct ? ffNLP::user_defined : ToLevel(autostructmode);
					if(A==bfgs || A==bfgs_affine_g) lh=ffNLP::do_nothing;
					_optim->BuildMatrixStructures(lh,lj,mmm);
					if(checkindex) _optim->EnableCheckStruct();
					
					SmartPtr<IpoptApplication> app = new IpoptApplication();
					
					//app->Options()->SetNumericValue("tol", 1e-10);
					if(nargs[9]) app->Options()->SetNumericValue("tol",GetAny<double>((*nargs[9])(stack)));
					if(nargs[10]) app->Options()->SetIntegerValue("max_iter",GetAny<long>((*nargs[10])(stack)));
					if(nargs[11]) app->Options()->SetNumericValue("max_cpu_time",GetAny<double>((*nargs[11])(stack)));
					//app->Options()->SetStringValue("hessian_approximation","limited-memory");
					if(A==bfgs || A==bfgs_affine_g) app->Options()->SetStringValue("hessian_approximation","limited-memory");
					app->Options()->SetStringValue("mu_strategy", "adaptive");
					app->Options()->SetStringValue("output_file", "ipopt.out");
					//app->Options()->SetStringValue("mehrotra_algorithm", "yes");
					
					ApplicationReturnStatus status;
					app->Initialize();
			
					// Ask Ipopt to solve the problem
					status = app->OptimizeTNLP(optim);
					
					if(lag_mul) *lag_mul = _optim->lambda_start;
					
					if (status == Solve_Succeeded) {
						printf("\n\n*** Ipopt succeeded \n");
					}
					else {
						printf("\n\n*** Ipopt failure!\n");
					}
					if(lag_mul) delete lag_mul;
					if(ffH) delete ffH;
					if(ffC) delete ffC;
					if(ffdC) delete ffdC;
					closetheparam.eval(stack); // clean memory 
					WhereStackOfPtr2Free(stack)->clean(); // FH mars 2005 
                                        lm.destroy(); // clean memory of LM 
					return cost; //SetAny<long>(0);  Modif FH  july 2005       
				}
				    
				operator aType () const { return atype<double>();} 
				
		};
		
		E_F0 * code(const basicAC_F0 & args) const {return new E_Ipopt(args,A,WC);}
		
		
		
		OptimIpopt(Case<no_assumption,with_constraints>) : 
			OneOperator(atype<double>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<KN<R> *>()),
			A(no_assumption),WC(with_constraints) {}
		OptimIpopt(Case<no_assumption,without_constraints>) :
			OneOperator(atype<double>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<KN<R> *>()),
			A(no_assumption),WC(without_constraints) {}
		OptimIpopt(Case<affine_g,with_constraints>) :
			OneOperator(atype<double>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Matrice_Creuse<R> *>(),atype<KN<R> *>()),
			A(affine_g),WC(with_constraints) {}
		OptimIpopt(Case<quadratic_f_affine_g,with_constraints>) :
			OneOperator(atype<double>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Matrice_Creuse<R> *>(),atype<Polymorphic*>(),atype<Matrice_Creuse<R> *>(),atype<KN<R> *>()),
			A(quadratic_f_affine_g),WC(with_constraints) {}
		OptimIpopt(Case<quadratic_f_affine_g,without_constraints>) :
			OneOperator(atype<double>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Matrice_Creuse<R> *>(),atype<KN<R> *>()),
			A(quadratic_f_affine_g),WC(without_constraints) {}
		OptimIpopt(Case<spurious_assumption,with_constraints>) :
			OneOperator(atype<double>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Matrice_Creuse<R> *>(),atype<Polymorphic*>(),atype<Polymorphic *>(),atype<KN<R> *>()),
			A(spurious_assumption),WC(with_constraints) {}
		OptimIpopt(Case<bfgs,with_constraints>) :
			OneOperator(atype<double>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<KN<R> *>()),
			A(bfgs),WC(with_constraints) {}
		OptimIpopt(Case<bfgs,without_constraints>) :
			OneOperator(atype<double>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<KN<R> *>()),A(bfgs),WC(without_constraints) {}
		OptimIpopt(Case<bfgs_affine_g,with_constraints>) :
			OneOperator(atype<double>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Polymorphic*>(),atype<Matrice_Creuse<R>*>(),atype<KN<R> *>()),
			A(bfgs_affine_g),WC(with_constraints) {}
			
};

void OptimIpopt::E_Ipopt::InitUNP()
{
	if(A==no_assumption && WC) {} //no unused named parameter
	if(A==no_assumption && !WC)					AddElements(unused_name_param,4,2,3,4,6);
	if(A==affine_g && WC)								AddElements(unused_name_param,1,4);
	if(A==quadratic_f_affine_g && WC)		AddElements(unused_name_param,4,4,5,7,8);
	if(A==quadratic_f_affine_g && !WC)	AddElements(unused_name_param,7,2,3,4,5,6,7,8);
	if(A==spurious_assumption)					AddElements(unused_name_param,12,0,1,2,3,4,5,6,7,8,9,10,11);
	if(A==bfgs && WC)										AddElements(unused_name_param,1,5);
	if(A==bfgs && !WC)									AddElements(unused_name_param,7,2,3,4,5,6,7,8);
	if(A==bfgs_affine_g)								AddElements(unused_name_param,4,4,5,7,8);
}


basicAC_F0::name_and_type  OptimIpopt::E_Ipopt::name_param[]= 
{
	{"lb",				&typeid(KN_<double>) },									//0
	{"ub",				&typeid(KN_<double>) },									//1
	{"clb",				&typeid(KN_<double>) },									//2
	{"cub",				&typeid(KN_<double>) },									//3
	{"structjac", &typeid(E_Array)},											//4
	{"structhess",&typeid(E_Array)},											//5
	{"lm",				&typeid(KN_<double>)},									//6
	{"autostruct",&typeid(long)},													//7
	{"checkindex",&typeid(bool)},													//8
	{"tol",				&typeid(double)},												//9
	{"maxiter",		&typeid(long)},													//10
	{"maxcputime",&typeid(double)}												//11
};



class Init { public:
  Init();
};

static Init init;

Init::Init()  
{
  Global.Add("IPOPT","(",new OptimIpopt(Case<no_assumption,with_constraints>()));
  Global.Add("IPOPT","(",new OptimIpopt(Case<no_assumption,without_constraints>())); 
	Global.Add("IPOPT","(",new OptimIpopt(Case<affine_g,with_constraints>()));
	Global.Add("IPOPT","(",new OptimIpopt(Case<quadratic_f_affine_g,with_constraints>()));
	Global.Add("IPOPT","(",new OptimIpopt(Case<quadratic_f_affine_g,without_constraints>()));
	Global.Add("IPOPT","(",new OptimIpopt(Case<spurious_assumption,with_constraints>()));
	Global.Add("IPOPT","(",new OptimIpopt(Case<bfgs,with_constraints>()));
	Global.Add("IPOPT","(",new OptimIpopt(Case<bfgs,without_constraints>()));
	Global.Add("IPOPT","(",new OptimIpopt(Case<bfgs_affine_g,with_constraints>()));
}









