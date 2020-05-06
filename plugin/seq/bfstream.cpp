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
#include <cstdio>
#include <unistd.h>

template< class T >
class Stream_b {
 public:
  Stream_b(T *ff) : f(ff) {}

  Stream_b(T **ff) : f(*ff) { ffassert(f); }

  Stream_b(const Stream_b &io) : f(io.f) {}

  T *f;
};

/*!
 * \brief pto_stream_b
 * \param f T**
 * \return Stream_b<T>
 */
template< class T >
Stream_b< T > pto_stream_b(T **f) {
  return Stream_b< T >(f);
}

/*!
 * \brief to_stream_b
 * \param f T*
 * \return Stream_b<T>
 */
template< class T >
Stream_b< T > to_stream_b(T *f) {
  return Stream_b< T >(f);
}

/*!
 * \brief Read
 * \param io Stream_b<istream> const &
 * \param data T *const &
 * \return istream *
 */
template< class T >
istream *Read(Stream_b< istream > const &io, T *const &data) {
  io.f->read(reinterpret_cast< char * >(data), sizeof(*data));
  return io.f;
}

/*!
 * \brief Read
 * \param io Stream_b<istream> const &
 * \param data KN<T> * const &
 * \return istream *
 */
template< class T >
istream *Read(Stream_b< istream > const &io, KN< T > *const &data) {
  long n;

  io.f->read(reinterpret_cast< char * >(&n), sizeof(long));
  if (verbosity > 0) cout << " read  n =" << n << " " << n * sizeof(sizeof(T)) << " " << endl;
  if (n != data->N( )) {
    data->resize(n);
  }

  T *p = *data;
  io.f->read(reinterpret_cast< char * >(p), n * sizeof(T));
  return io.f;
}

/*!
 * \brief Write
 * \param io Stream_b<ostream> const &
 * \param data KN<T> * const
 * \return ostream *
 */
template< class T >
ostream *Write(Stream_b< ostream > const &io, KN< T > *const &data) {
  T *p = *data;
  long n = data->N( );

  if (verbosity > 0) cout << " write n =" << n << " " << n * sizeof(T) << " " << p << endl;
  io.f->write(reinterpret_cast< const char * >(&n), sizeof(long));
  io.f->write(reinterpret_cast< const char * >(p), n * sizeof(T));
  return io.f;
}

/*!
 * \brief Write
 * \param io Stream_b<ostream> const &
 * \param data T * const &
 * \return ostream *
 */
template< class T >
ostream *Write(Stream_b< ostream > const &io, T *const &data) {
  io.f->write(reinterpret_cast< const char * >(data), sizeof(*data));
  return io.f;
}

/*!
 * \brief Write
 * \param io Stream_b<ostream> const &
 * \param data T * const &
 * \return ostream *
 */
template< class T >
ostream *Write(Stream_b< ostream > const &io, T const &data) {
  io.f->write(reinterpret_cast< const char * >(&data), sizeof(data));
  return io.f;
}

/*!
 * \brief init k
 * \param io Stream_b<ostream> const &
 * \param data T * const &
 * \return ostream *
 */
template< class K >
void initK( ) {
  typedef Stream_b< ostream > OB;
  typedef Stream_b< istream > IB;
  Add< IB >("(", "", new OneOperator2_< istream *, IB, K * >(Read));
  Add< OB >("(", "", new OneOperator2_< ostream *, OB, K * >(10,Write));
  Add< OB >("(", "", new OneOperator2_< ostream *, OB, K >(Write));
  Add< IB >("(", "", new OneOperator2_< istream *, IB, KN< K > * >(Read));
  Add< OB >("(", "", new OneOperator2_< ostream *, OB, KN< K > * >(Write));
}

static void inittt( ) {
  typedef Stream_b< ostream > OB;
  typedef Stream_b< istream > IB;
  Dcl_Type< OB >( );
  Dcl_Type< IB >( );

  Add< istream ** >("read", ".", new OneOperator1< IB, istream ** >(pto_stream_b< istream >));
  Add< ostream ** >("write", ".", new OneOperator1< OB, ostream ** >(pto_stream_b< ostream >));
  initK< long >( );
  initK< double >( );
  initK< complex< double > >( );
}

LOADFUNC(inittt);
