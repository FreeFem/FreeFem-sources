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
#include <cstdint>
template<size_t sz>
struct swap_bytes
{
    static   char *  swap(char *t)
    {
        throw std::out_of_range("data size");
    }
};
template<>
struct swap_bytes<1>
{
    static inline char *  swap(char *t)
    {
        return t;
    }
};
template<>
struct swap_bytes<2>
{
    static inline char *  swap(char *t)
    {
        std::swap(t[0],t[1]);
        return t;
    }
};
void dumpb(char *t,int sz)
{
    for(int k=0; k< sz; ++k)
    cout << (int) t[k] << " ";
    cout << "\n";
}

template<>
struct swap_bytes<4>
{
    static inline char *  swap(char *t)
    {
        dumpb(t,4);
        std::swap(t[0],t[3]);
        std::swap(t[1],t[2]);
        dumpb(t,4);
        return t;
    }
};
template<>
struct swap_bytes<8>
{
    static   inline char *  swap(char *t)
    {
        std::swap(t[0],t[7]);
        std::swap(t[1],t[6]);
        std::swap(t[2],t[5]);
        std::swap(t[3],t[4]);
        return t;
    }
};

bool islittleendian() { uint8_t i=1,ii = (i <<8); return  ii==256;}
bool isbigendian() {  uint8_t i=1,ii = (i <<8);return  ii!=256;}

template<typename T> T swapbyte(T val) { swap_bytes<sizeof(T)>::swap( reinterpret_cast< char * >( &val));return val;}

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

template< class T,class TW >
ostream *Write( ostream * f, T data) {
    TW dataw=data;
    f->write(reinterpret_cast< const char * >(&dataw), sizeof(T));
    return f;
}

template< class T >
istream * Reada( istream * f, KN_<T>  a) {
    if( a.contiguous())
        f->read(reinterpret_cast<  char * >((T*) a), sizeof(T)*a.N());
    else
        for(int i=0; i<a.N();++i)
    f->read(reinterpret_cast<  char * >(& a[i] ), sizeof(T));
    return f;
}

template< class T, class TR >
istream * Read( istream * f, T * data) {
    TR datar;
    f->read(reinterpret_cast<  char * >(&datar), sizeof(TR));
    *data = datar;
    return f;
}
template< class T, class TR >
istream * readswapbyte( istream * f, T * data) {
    TR datar;
    f->read(reinterpret_cast<  char * >(&datar), sizeof(TR));
    datar=  swapbyte(datar);
    *data = datar;
    return f;
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
    
    Global.Add("read", "(",
               new OneOperator2<istream *,istream *,Complex*>(Read<Complex,Complex>) ,
               new OneOperator2<istream *,istream *,double*>(Read<double,double>) ,
               new OneOperator2<istream *,istream *,long*>(Read<long,long>),
               new OneOperator2<istream *,istream *,KN_<Complex> >(Reada<Complex>) ,
               new OneOperator2<istream *,istream *,KN_<double> >(Reada<double>) ,
               new OneOperator2<istream *,istream *,KN_<long> >(Reada<long>)
               );
    
    Global.Add("readswapbyte", "(",
               new OneOperator2<istream *,istream *,double*>(readswapbyte<double,double>) ,
               new OneOperator2<istream *,istream *,long*>(readswapbyte<long,long>)
               );
    Global.Add("readint", "(",
               new OneOperator2<istream *,istream *,long*>(Read<long,int>)
               );
    Global.Add("readfloat", "(",
               new OneOperator2<istream *,istream *,double*>(Read<double,float>)
               );
    Global.Add("writeint", "(",
               new OneOperator2<ostream *,ostream *,long>(Write<long,int>)
               );
    Global.Add("writefloat", "(",
               new OneOperator2<ostream *,ostream *,double>(Write<double,float>)
               );
    Global.Add("readintswapbyte", "(",
               new OneOperator2<istream *,istream *,long*>(readswapbyte<long,int>)
               );
    
    if(verbosity> 9)
        cout << "\n islittleendian =" <<  islittleendian() << endl
        << " or  isbigeendian =" <<  isbigendian() << endl;
}

LOADFUNC(inittt);
