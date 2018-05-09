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

template<class T> struct  Resize {T *v;
	                              Resize (T *vv): v(vv) {}
};

template<class T> T*resize1 (const Resize<T> &t, const long &n) {
	t.v->resize(n);
	return t.v;
}

template<class T> T*resizeandclean1 (const Resize<T> &t, const long &n) {
	int nn = t.v->N();	// old size

	for (int i = n; i < nn; i++) delete (*t.v)[i];	// clean

	t.v->resize(n);

	for (int i = nn; i < n; i++) (*t.v)[i] = 0;

	return t.v;
}

template<class T> T*resize2 (const Resize<T> &t, const long &n, const long &m) {
	t.v->resize(n, m);
	return t.v;
}

template<class T> Resize<T> to_Resize (T *v) {return Resize<T>(v);}

template<class T> struct  Resize1 {T v;
	                               Resize1 (T vv): v(vv) {}
};
template<class T> Resize1<T> to_Resize1 (T v) {return Resize1<T>(v);}

