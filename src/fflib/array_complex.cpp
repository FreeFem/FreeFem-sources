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
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

#include "array_tlp.hpp"
#include "array_init.hpp"

/*!
 * \brief Get_C2R
 * Not used ?
 */
// KN_<double> Get_C2R(KN_<complex<double> > vc) {
//   ffassert(vc.step == 1);
//   complex<double> * pc = vc; // table pointer
//   double *pr = static_cast<double*>(static_cast<void*>(pc));
//   return KN_<double>(pr, vc.N()*2);
// }

/*!
 * \brief Separate real and imaginary part of complex array
 * \param vc Complex array, size n
 * \return Real or Imaginary part of the complex array, depending on offset (0: real, 1: imaginary)
 */
template<int offset>
KN_<double> Get_C2R_(KN_<complex<double> > vc) {
  complex<double> * pc = vc; // table pointer
  double *pr = static_cast<double*>(static_cast<void*>(pc));
  return KN_<double>(pr+offset, vc.N(), vc.step*2);
}

/*!
 * \brief Separate real and imaginary part of complex array
 * \param vc Complex array, size n*m
 * \return Real or Imaginary part of the complex array, depending on offset (0: real, 1: imaginary)
 */
template<int offset>
KNM_<double> Get_C2R_(KNM_<complex<double> > vc) {
  complex<double> * pc = vc; // table pointer
  double *pr = static_cast<double*>(static_cast<void*>(pc));
  return KNM_<double>(pr+offset, ShapeOfArray(vc.N(), vc.step*2), vc.shapei, vc.shapej);
}

/*!
 * \brief Initialization
 */
KN_<Complex> rmeps(KN_rmeps<Complex> p,double eps)
{
    for(int i=0;i<p.v.N();++i)
    {   int c=0;double a =p.v[i].real(), b=p.v[i].imag();
        if(abs(a)<eps) c=1,a=0.;
        if(abs(b)<eps) c=1,b=0.;
        if(c) p.v[i]=Complex(a,b);
    }
    return p.v;}//
void initArrayDCLComplex() {
  ArrayDCL<Complex>();
  Dcl_Type<Mul_KNMh_KN_<Complex> >();
}

Complex square(const Complex &x){ return x*x; }
// double imagc(const Complex &x){ return imag(x); }  // Commented out in interface
// double realc(const Complex &x){ return realc(x); }

void initArrayOperatorComplex() {
  typedef Complex K;
  typedef const K & KK;
     Dcl_Type< KN_rmeps<K> > ();
  typedef double R;
  typedef const R & RR;
  ArrayOperator<Complex, long>();
  ArrayOperatorF<Complex, const Complex&>();
  Add<KN<K> *>("rmeps",".",new OneOperator1<KN_rmeps<K>,KN_<K> >(build_rmeps));
  Add<KN_rmeps<K> >("(","",new OneOperator2<KN_<K>,KN_rmeps<K>,double>(rmeps));

  // Take the real or imaginary part of a complex array
  Add<KN<Complex> *>("im", ".", new OneOperator1<KN_<double>, KN_<Complex> >(Get_C2R_<1>, atype<KN<Complex>* >()));
  Add<KN_<Complex> >("im", ".", new OneOperator1<KN_<double>, KN_<Complex> >(Get_C2R_<1>, atype<KN_<Complex> >()));
  Add<KN<Complex> *>("re", ".", new OneOperator1<KN_<double>, KN_<Complex> >(Get_C2R_<0>, atype<KN<Complex>* >()));
  Add<KN_<Complex> >("re", ".", new OneOperator1<KN_<double>, KN_<Complex> >(Get_C2R_<0>, atype<KN_<Complex> >()));

  Add<KNM<Complex> *>("im", ".", new OneOperator1<KNM_<double>, KNM_<Complex> >(Get_C2R_<1>, atype<KNM<Complex>* >()));
  Add<KNM_<Complex> >("im", ".", new OneOperator1<KNM_<double>, KNM_<Complex> >(Get_C2R_<1>, atype<KNM_<Complex> >()));
  Add<KNM<Complex> *>("re", ".", new OneOperator1<KNM_<double>, KNM_<Complex> >(Get_C2R_<0>, atype<KNM<Complex>* >()));
  Add<KNM_<Complex> >("re", ".", new OneOperator1<KNM_<double>, KNM_<Complex> >(Get_C2R_<0>, atype<KNM_<Complex> >()));

  Global.Add("square", "(", new OneOperator1F_KN_<F_KN_<K, K, K, KK>, K, KK, KN_<K> >(square));
  Global.Add("conj", "(", new OneOperator1F_KN_<F_KN_<K, K, K, KK>,K, KK, KN_<K> >(conj));
  // Global.Add("imag", "(", new OneOperator1F_KN_<F_KN_<R, R, K, KK>, R, KK, KN_<K> >(imagc));
  // Add<KN_<double> >("<--", "(", new OneOperator1F_KN_<F_KN_<R, R, K, KK>, R, KK, KN_<K> >(realc));
  TheOperators->Add("*", new OneBinaryOperator<Op2_2p_<Mul_KNMh_KN_<K>, Transpose<KNM<K>*>, KN<K>*> >); // A'*b
  TheOperators->Add("=", new OneBinaryOperator<init_eqarray<KN<K>, Mul_KNMh_KN_<K> > >);
  TheOperators->Add("<-", new OneBinaryOperator<init_eqarray<KN<K>, Mul_KNMh_KN_<K> > >);
  map_type[typeid(SetArray<K>).name()]->AddCast(
    new E_F1_funcT<SetArray<K>, SetArray<long> >(Cast<SetArray<K>, SetArray<long> >),
    new E_F1_funcT<SetArray<K>, SetArray<double> >(Cast<SetArray<K>, SetArray<double> >)
  );
  Global.Add("toCarray", "(", new OneOperator_2KN_<Complex>);
}
