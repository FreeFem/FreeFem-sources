// to compile ff-c++ pipe.cpp
//  warning do not compile under windows...
#include "ff++.hpp"
#include <ext/stdio_filebuf.h>
#include <cstdio> 
#include <unistd.h>
//#include "pstream.h"

typedef __gnu_cxx::stdio_filebuf<char> stdiofilebuf;
long ffsleep(long  s) { return sleep(s);}
long ffusleep(long  s) { return usleep(s);}
struct pstream {
    FILE * f; 
    stdiofilebuf * fb;
    ostream *os;
    istream *is;
    pstream(FILE *ff,std::ios_base::openmode mode)
     : f(ff),fb(new stdiofilebuf(f,mode) ) ,os(0),is(0)
    {  if(verbosity>10) cout << " mode " << mode << endl;
        if(mode & ios_base::in )  is = new istream(fb);
        if(mode & ios_base::out ) os = new ostream(fb);
        if(verbosity>10)  cout << is << " " << os << " ******* " << endl;
    }
    ~pstream() {
        if(f) pclose(f);
        if(os) delete os;
        if(is) delete is;
        if(fb) delete(fb);
        f=0;
        os=0;
        is=0;
        fb=0;
}
};

//typedef redi::pstream pstream;
typedef std::string string;

static pstream **  pstream_init(pstream **const & p,string * const & a,string * const & b)
{
    string mode= b ? *b: "w";
    if (mode.length() ==0) mode = "wr";
    std::ios_base::openmode  om = ios_base::in| ios_base::out;
    if (mode == "r+" )      om =ios_base::in| ios_base::out;
    else if( mode == "w"  ) om = ios_base::out;
    else if(  mode == "r" ) om =  ios_base::in;
    else  ExecError("Invalide mode pstream r,r+,w ");
   if(verbosity>10)  *ffapi::cout()  << "pstream_init: om " << om << "(" <<ios_base::in << ios_base::out << ") mode:"
                    << mode << " '" << *a <<"'"<<  endl;
#ifdef WIN32
    FILE * pp =_popen(a->c_str(),mode.c_str());
#else    
    FILE * pp =popen(a->c_str(),mode.c_str());
#endif
    *p =  new pstream(pp,om);

    if ( !*p || !pp) {
        cerr << " Error openning pipe  " << *a << endl;
        ExecError("Error openning pipe");}
   
    return p;
};

static pstream **  pstream_init(pstream **const & p,string * const & a) 
{
    return pstream_init(p,a,0);
}; 


AnyType pstream2o(Stack,const AnyType &a) {
    pstream* p = *PGetAny<pstream*>(a);
    ffassert(p->os);
    return   SetAny<ostream*>(p->os);}
    
AnyType pstream2i(Stack,const AnyType &a) {
        pstream* p = *PGetAny<pstream*>(a);
        ffassert(p->is);
        return   SetAny<istream*>(p->is);}

class istream_good { public:
    istream_good(istream * ff) :f(ff) {}
    istream * f;
    operator bool () const {return f->good();}
};
long cflush(pstream ** ppf) 
{
  pstream & f = **ppf;
  if( f.os ) f.os->flush();
  if( f.f) fflush(f.f);
  return  0; 
}; 
inline istream_good to_istream_good(pstream **f){ ffassert((**f).is) ; return istream_good((**f).is);}
inline bool get_eof(pstream ** p){ return (**p).is ? (**p).is->eof(): EOF;}
void inittt()
{
  Dcl_TypeandPtr<pstream*>(0,0,::InitializePtr<pstream*>,::DeletePtr<pstream*>);
  atype<istream* >()->AddCast( new E_F1_funcT<istream*,pstream**>(pstream2i));
  atype<ostream* >()->AddCast( new E_F1_funcT<ostream*,pstream**>(pstream2o));
  TheOperators->Add("<-",new OneOperator2_<pstream**,pstream**,string*>(pstream_init) );
  TheOperators->Add("<-",new OneOperator3_<pstream**,pstream**,string*,string*>(pstream_init) );
  zzzfff->Add("pstream",atype< pstream ** >());
  Add<pstream**>("good",".",new OneOperator1<istream_good,pstream**>(to_istream_good));
  Add<pstream**>("eof",".",new OneOperator1<bool,pstream**>(get_eof));
  Global.Add("flush","(",new OneOperator1<long,pstream **> ( cflush)) ; 
  Global.Add("sleep","(",new OneOperator1<long,long> ( ffsleep)) ; 
  Global.Add("usleep","(",new OneOperator1<long,long> ( ffusleep)) ; 
#ifdef WIN32
   Global.New("onWIN32",CConstant<bool>(true));
#else    
  Global.New("onWIN32",CConstant<bool>(false));
#endif
}

LOADFUNC(inittt);
