// to compile ff-c++ bstream.cpp
//  warning do not compile under windows...
#include "ff++.hpp"
#include <cstdio>
#include <unistd.h>

template<class T>
class Stream_b { public:
    Stream_b(T * ff) :f(ff) {}
    Stream_b(T ** ff) :f(*ff) {ffassert(f); }
    Stream_b(const Stream_b &io): f(io.f) {}
    
    T * f;
};

template<class T>
 Stream_b<T>  pto_stream_b(T **f){ return Stream_b<T>(f);}
template<class T>
Stream_b<T>  to_stream_b(T *f){ return Stream_b<T>(f);}

template<class T>
istream * Read(Stream_b<istream> const &   io, T * const & data ) {
    io.f->read(reinterpret_cast<char *>(data),sizeof(*data));
    return io.f; }

template<class T>
istream * Read(Stream_b<istream> const &   io, KN<T> * const & data ) {
     long n;
    io.f->read(reinterpret_cast<char *>(&n),sizeof(long));
    cout << " read  n =" << n << " " << n*sizeof(sizeof(T)) << " "  <<  endl;
    if( n != data->N()) data->resize(n);
    T*  p = *data;
    io.f->read(reinterpret_cast<char *>(p),n*sizeof(T));
    return io.f; }
template<class T>
ostream * Write(Stream_b<ostream> const &   io, KN<T> * const & data ) {
    T*  p = *data;
    long n=data->N();
    cout << " write n =" << n << " " << n*sizeof(T) << " " << p <<  endl;
    io.f->write(reinterpret_cast<const char *>(&n),sizeof(long));
    io.f->write(reinterpret_cast<const char *>(p),n*sizeof(T));
    return io.f; }

template<class T>
ostream * Write(Stream_b<ostream> const &   io, T * const & data ) {
    io.f->write(reinterpret_cast<const char *>(data),sizeof(*data));
    return io.f; }
template<class T>
ostream * Write(Stream_b<ostream> const &   io, T  const & data ) {
    io.f->write(reinterpret_cast<const char *>(&data),sizeof(data));
    return io.f; }
template <class K>
void initK()
{
    typedef Stream_b<ostream> OB;
    typedef Stream_b<istream> IB;
    Add<IB>("(","",new OneOperator2_<istream *,IB,K *>(Read));
    Add<OB>("(","",new OneOperator2_<ostream *,OB,K *>(Write));
    Add<OB>("(","",new OneOperator2_<ostream *,OB,K >(Write));
    Add<IB>("(","",new OneOperator2_<istream *,IB,KN<K> *>(Read));
    Add<OB>("(","",new OneOperator2_<ostream *,OB,KN<K> * >(Write));
    
}
static void inittt()
{
    typedef Stream_b<ostream> OB;
    typedef Stream_b<istream> IB;
    Dcl_Type< OB>  ();
    Dcl_Type< IB> ();
    
    Add<istream**>("read",".",new OneOperator1<IB,istream**>(pto_stream_b<istream>));
    Add<ostream**>("write",".",new OneOperator1<OB,ostream**>(pto_stream_b<ostream>));
    initK<long>();
    initK<double>();
    initK<complex<double> >();
    
/*
     Add<IB>("(","",new OneOperator2_<istream *,IB,double *>(Read));
    Add<OB>("(","",new OneOperator2_<ostream *,OB,double *>(Write));
    Add<OB>("(","",new OneOperator2_<ostream *,OB,double >(Write));
 */
    
}

LOADFUNC(inittt);
