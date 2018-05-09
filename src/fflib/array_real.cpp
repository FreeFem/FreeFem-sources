/*
 * This file is part of FreeFem++.
 *
 * FreeFem++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FreeFem++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

// -*- Mode : c++ -*-
//
// SUMMARY  :
// USAGE    :
// ORG      :
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

#include "array_tlp.hpp"
#include "array_init.hpp"

double square (double x) {return x * x;}

void initArrayDCLdouble () {
// ArrayOperator<long>();
	ArrayDCL<double>();
}

// template<class A,class B>  A Build(B b) {  return A(b);}
void initArrayOperatordouble () {
	ArrayOperator<double, long>();
	ArrayOperatorF<double, double>();
	typedef double K;
	typedef double KK;
	Dcl_Type<QuantileKN<K>>();
	Global.Add("abs", "(", new OneOperator1F_KN_<F_KN_<K, K, K, KK>, K, KK, KN_<K>>(fabs));
	Global.Add("acos", "(", new OneOperator1F_KN_<F_KN_<K, K, K, KK>, K, KK, KN_<K>>(acos));
	Global.Add("asin", "(", new OneOperator1F_KN_<F_KN_<K, K, K, KK>, K, KK, KN_<K>>(asin));
	Global.Add("atan", "(", new OneOperator1F_KN_<F_KN_<K, K, K, KK>, K, KK, KN_<K>>(atan));
	Global.Add("floor", "(", new OneOperator1F_KN_<F_KN_<K, K, K, KK>, K, KK, KN_<K>>(floor));
	Global.Add("ceil", "(", new OneOperator1F_KN_<F_KN_<K, K, K, KK>, K, KK, KN_<K>>(ceil));

	Global.Add("square", "(", new OneOperator1F_KN_<F_KN_<K, K, K, KK>, K, KK, KN_<K>>(square));
	Add<KN_<double>>("sort", ".", new OneOperator1_<KN_<K>, KN_<K>>(SortKn<K, KN_<K>> ));
	// Add<KN<double> >("sort",".",new OneOperator1_<KN<K>,KN<K> >(SortKn<K, KN<K> >));
	Add<KN<double> *>("sort", ".", new OneOperator1_<KN<K> *, KN<K> *>(SortpKn<K> ));

	// Add<KN_<double> >("sort",".",new OneOperator2_<KN_<K>,KN_<K>,KN_<long> >(SortKn<K,long, KN_<K>,KN_<long>  >));
	// Add<KN<double> >("sort",".",new OneOperator2_<KN<K>,KN<K>,KN_<long> >(SortKn<K,long, KN<K>,KN_<long> >));
	Global.Add("sort", "(", new OneOperator2_<KN<K> *, KN<K> *, KN<long> *>(SortpKn2<K, long> ));

	Add<KN_<double>>("quantile", ".", new OneOperator1<QuantileKN<K>, KN_<K>>(Build<QuantileKN<K>, KN_<K>> ));
	// Add<KN<double> >("quantile",".",new OneOperator1<QuantileKN<K>,KN_<K> >(Build<QuantileKN<K>,KN_<K> >, atype<KN<K> >()));
	// Add<KN_<double> * >("quantile",".",new OneOperator1<QuantileKN<K>,KN_<K> >(Build<QuantileKN<K>,KN_<K> >, atype<KN<K> >()));
	Add<KN<double> *>("quantile", ".", new OneOperator1<QuantileKN<K>, KN_<K>>(Build<QuantileKN<K>, KN_<K>>, atype<KN<K> *>()));
	Add<QuantileKN<K>>("(", "", new OneOperator2_<K, QuantileKN<K>, double>(Quantile<K> ));

	map_type[typeid(SetArray<K>).name()]->AddCast(
		new E_F1_funcT<SetArray<K>, SetArray<long>>(Cast<SetArray<K>, SetArray<long>> )

		);
	Global.Add("toRarray", "(", new OneOperator_2KN_<double> );
	Add<KN<K> *>("imin", ".", new OneOperator1<long, KN<K> *>(get_imin));
	Add<KN<K> *>("imax", ".", new OneOperator1<long, KN<K> *>(get_imax));
	Add<KNM<K> *>("imin", ".", new OneOperator1<long, KNM<K> *>(get_imin));	// Add april 2018 FH
	Add<KNM<K> *>("imax", ".", new OneOperator1<long, KNM<K> *>(get_imax));	// Add april 2018 FH
	Add<KNM<K> *>("jmin", ".", new OneOperator1<long, KNM<K> *>(get_jmin));	// Add april 2018 FH
	Add<KNM<K> *>("jmax", ".", new OneOperator1<long, KNM<K> *>(get_jmax));	// Add april 2018 FH
	Global.Add("ijmax", "(", new OneOperator3_<NothingType, KNM<K> *, long *, long *>(get_ijmax));	// Add april 2018 FH
	Global.Add("ijmin", "(", new OneOperator3_<NothingType, KNM<K> *, long *, long *>(get_ijmin));	// Add april 2018 FH
	// madd FH. march 2015 ...

	Global.Add("Unique", "(", new Unique<K, K> );
	Global.Add("Unique", "(", new Unique<K, long> );

// ArrayDCL<long>();
}

