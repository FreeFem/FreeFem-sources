/*
 *  PlotStream.hpp
 *
 *  Created by Frederic Hecht on 27/10/08.
 *  Copyright 2008 UPMC.
 *
 */

#include <cstdio>
#ifdef _WIN32
#include <fcntl.h>
#endif
#include "endian.hpp"

#ifndef PLOT_STREAM_HPP_
#define PLOT_STREAM_HPP_

//FFCS visualization stream redirection
#include "ffapi.hpp"

using  Fem2D::Mesh;
using  Fem2D::Mesh3;
namespace Fem2D {
    
}


class PlotStream 
{
public:
    
  FILE * TheStream;

  PlotStream(FILE *Stream) :TheStream(Stream) { }
  operator bool() const { return TheStream;}
  // datatype mush be < 0 to have no collistion with arg number. 
  // FFCS: <<PlotStream::datatype>>
    enum datatype { dt_meshes=-1,dt_plots=-2,dt_endplot=-3,dt_meshes3=-10,dt_endarg=99999,dt_newplot=-5, dt_meshesS=-12, dt_meshesL=-13};    // modif order
  
  //  dt_plots3=-11,dt_plotsS=-13,dt_plots3S=-14
  //FFCS:need to send control data to FFCS at least once per plot
  void SendNewPlot() {  ffapi::newplot();write((long )dt_newplot); set_binary_mode(); }
  void SendEndArgPlot() {write((long )dt_endarg); }
  //FFCS:redirect visualization stream
  void SendEndPlot() { write((long )dt_endplot);ffapi::ff_fflush(TheStream); set_text_mode() ;}
  void SendPlots() { write((long )dt_plots); }
  void SendMeshes() { write((long )dt_meshes);}
  void SendMeshes3() { write((long )dt_meshes3);}
  void SendMeshesS() { write((long )dt_meshesS);}
  void SendMeshesL() { write((long )dt_meshesL);}
  //FFCS: divert stream to FFCS
  void write(const void *data,size_t l) {ffapi::ff_fwrite(data,1,l,TheStream);}

  PlotStream& write(const bool& bb) {bool b=w_endian(bb);write(reinterpret_cast<const void *> (&b),sizeof(bool));return *this;}
  PlotStream& write(const long long& bb) {long long  b=w_endian(bb);write(reinterpret_cast<const void *> (&b),sizeof(long long));return *this;}
  PlotStream& write(const long& bb) { // always write 8 bits for  a long FH. 
    long long ll=bb;ll=w_endian(ll);write(reinterpret_cast<const void *> (&ll),sizeof(long long));
    return *this;}
  PlotStream& write(const int& bb) {int b=w_endian(bb);write(reinterpret_cast<const void *> (&b),sizeof(int));return *this;}
  PlotStream& write(const double& bb) {double b=w_endian(bb);write(reinterpret_cast<const void *> (&b),sizeof(double));return *this;}
  PlotStream& write(const complex<double>& bb) {return write(bb.real()),write(bb.imag());}
  PlotStream &write(const Fem2D::R1 & P) { return write(P.x);}
  PlotStream &write(const Fem2D::R2 & P) { return write(P.x),write(P.y);}
  PlotStream &write(const Fem2D::R3 & P) { return write(P.x),write(P.y),write(P.z);}

  PlotStream& write(const string& b) {  
    int l=b.size();
    write(l);
    write(b.data(),l);
    return *this;
  }
  
  void set_text_mode() 
  {    
  //FFCS:visualization stream redirection
    ffapi::wintextmode(TheStream);
  }
  void set_binary_mode()
  {    
  //FFCS:visualization stream redirection
    ffapi::winbinmode(TheStream);
  }
 
  PlotStream & operator << (const bool& b)    { return write(b); }
  PlotStream & operator << (const long& b)    { return write(b); }        
  PlotStream & operator << (const long long & b)   { return write(b); }        
  PlotStream & operator << (const int& b)      { return write(b); }        
  PlotStream & operator << (const double& b)   { return write(b); }
  PlotStream & operator << (const complex<double>& b)   { return write(b); }
  PlotStream & operator << (const string& s)   { return write(s); }
  PlotStream & operator << (const string* s)   { ffassert(s); return write(*s); }
    template<class T>    
    PlotStream & operator << (const KN_<T>& b)
    {
	long n=b.N();
	write(n);
	// cout << "PlotStream :<<  n " << n << endl;

	for (int i=0;i<n;++i)
	    write(b[i]);
	return *this;
    }


  PlotStream & operator << (const Mesh& Th) {
      /*
    Serialize s=Th.serialize();
    long  n=s.size();
    write( n );
    write(s,s.size());*/
    GSave2(TheStream , Th) ; 
    return *this;}

  PlotStream & operator << (const Fem2D::Mesh2& Th) { Th.GSave(TheStream); return *this;}
  PlotStream & operator << (const Fem2D::Mesh3& Th) { Th.GSave(TheStream); return *this;}
  PlotStream & operator << (const Fem2D::MeshS& Th) { Th.GSave(TheStream); return *this;}
  PlotStream & operator << (const Fem2D::MeshL& Th) { Th.GSave(TheStream); return *this;}
    
  void read( void *data,size_t l) {
	char * p= (char*)data;
	for(int i=0;i<l;++i)	
	  *p++ = (char) getc(TheStream);
	//	fread(  p++,1,1,TheStream);
	//       read(data,l);
}
  //FFCS:visualization stream redirection
  bool good() const {return ffapi::ff_ferror(TheStream)==0;}
  void GetNewPlot() { get(dt_newplot) ; set_binary_mode();}
  void GetEndArgPlot() {get(dt_endarg); }
  void GetEndPlot() {get(dt_endplot); set_text_mode();}
  void GetPlots() { get(dt_plots); }
  void GetMeshes() { get(dt_meshes);}
 //void GetMeshesS() { get(dt_meshesS);}
 /* bool GetMeshes3() { long tt; read(tt);
      if(tt== dt_meshes3) return true;
      else if (tt== dt_plots) return false;
      cout << " Error Check :  get " << tt << " == wait for  "<< dt_meshes3 << " or "<< dt_plots << endl;
      ffassert(0);
     }*/
    
   inline int GetMeshes3() { long tt; read(tt);
      
       if(tt== dt_meshes3) return 0;
        else if (tt== dt_meshesS) return 1;
       else if (tt== dt_meshesL) return 2;
        else if (tt== dt_plots) return 3;    /// modif dans ff gl
     
        cout << " Error Check :  get " << tt << " == wait for  "<< dt_meshes3 << " or "<< dt_meshesS << " or "<< dt_meshesL << " or " << dt_plots << endl;
        ffassert(0);
    }

  void get(datatype t) { long tt; read(tt);
    if( tt !=(long) t) 
      cout << " Error Check :  get " << tt << " == wait for  "<< t << endl; 
    ffassert(tt==(long) t);}
  //FFCS:visualization stream redirection
    bool eof() {return ffapi::ff_feof(TheStream);}
  PlotStream& read( bool& b) {read(reinterpret_cast< void *> (&b),sizeof(bool));  b=r_endian(b);return *this;}
  PlotStream& read( long long& b) {read(reinterpret_cast< void *> (&b),sizeof(long long)); b=r_endian(b);return *this;}
  PlotStream& read( long& b) { long long l;
    read(reinterpret_cast< void *> (&l),sizeof(long long));
    l=r_endian(l);
    b=(long) l;
    if(( b-(long) l) !=0)
      { cout << " err err read long : error " << b << " !=  " << l << endl;
	assert( (b-(long) l)==0);}
    return *this;}
  PlotStream& read( int& b) {read(reinterpret_cast< void *> (&b),sizeof(int)); b=r_endian(b);return *this;}
  PlotStream& read( double& b) {read(reinterpret_cast< void *> (&b),sizeof(double)); b=r_endian(b);return *this;}
  PlotStream &read( complex<double> &C) { double  re,im; read(re); read(im); C=complex<double>(re,im); return 
 *this;}
  PlotStream &read( Fem2D::R1 & P) { return read(P.x);}
  PlotStream &read( Fem2D::R2 & P) { return read(P.x),read(P.y);}
  PlotStream &read( Fem2D::R3 & P) { return read(P.x),read(P.y),read(P.z);}
  PlotStream& read( string& b) {  	
    int l;
    read(l);
    b.resize(l);
    read(& (b[0]),l);
    return *this;
  }
  
  
  
  PlotStream & operator >> ( bool& b)      { return read(b); }
  PlotStream & operator >> ( long& b)      { return read(b); }        
  PlotStream & operator >> ( long long& b) { return read(b); }        
  PlotStream & operator >> ( int& b)       { return read(b); }        
  PlotStream & operator >> ( double& b)    { return read(b); }
  PlotStream & operator >> ( complex<double>& b)    { return read(b); }
  PlotStream & operator >> ( string& s)    { return read(s); }
  PlotStream & operator >> ( string *& s) 
  { if(!s) s= new string(); return read(*s);
    // cout << " fread string " << s <<endl;
  }
  
  PlotStream & operator >> ( Mesh *& Th)
  {
    long n;
    read(n);
    Serialize s(n,Mesh::magicmesh);
    read(s,n );
    Th= new Mesh(s);
    return *this;
  }
 
 
  template<class T>  
  PlotStream & operator >> ( KN<T>& b)
  {
    long n;
    read(n);
   // cout << "PlotStream >>  : n " << n << endl;
      //  read empty array .... (n ==0)
    if( ! b.N() && n) b.init(n);//   if n ==0 do nothing ..  Add FH nov. 2016
    ffassert( b.N()==n); 
    for (int i=0;i<n;++i)
      read(b[i]);
    return *this;
  }
  //  PlotStream & operator << (const Mesh3& Th);   
  //  PlotStream & operator >> ( Mesh3 *& Th);
    PlotStream &  operator >> ( Fem2D::Mesh3 *& Th) {	Th= new Fem2D::Mesh3(TheStream); return *this;}
    PlotStream &  operator >> ( Fem2D::Mesh2 *& Th) {	Th= new Fem2D::Mesh2(TheStream); return *this;}
    PlotStream &  operator >> ( Fem2D::MeshS *& Th) {    Th= new Fem2D::MeshS(TheStream); return *this;}
    PlotStream &  operator >> ( Fem2D::MeshL *& Th) {    Th= new Fem2D::MeshL(TheStream); return *this;}
    
    // ---   I also write the type .. to skip data if we  need  to skip data 
    // just change   >> and <<  by :  <= and >= 
    PlotStream & operator <= (const bool& b)    { return write(1),write(b); }
    PlotStream & operator <= (const long& b)    { return write(2),write(b); }        
    PlotStream & operator <= (const long long & b)   { return write(3),write(b); }        
    PlotStream & operator <= (const int& b)      { return write(4),write(b); }        
    PlotStream & operator <= (const double& b)   { return write(5),write(b); }
    PlotStream & operator <= (const string& s)   { return write(6),write(s); }
    PlotStream & operator <= (const string* s)   { return write(6),write(*s); }
    template<class T>    
    PlotStream & operator <= (const KN_<T>& b)   { return write(10),write((int) sizeof(T)),operator<<(b);}
    
    PlotStream & operator >= ( bool& b)       { return readc(1)>>b; }
    PlotStream & operator >= ( long& b)       { return readc(2)>>b; }        
    PlotStream & operator >= ( long long & b) {return  readc(3)>>b; }        
    PlotStream & operator >= ( int& b)        { return readc(4)>>b; }        
    PlotStream & operator >= ( double& b)     { return readc(5)>>b; }
    PlotStream & operator >= ( string& s)     { return readc(6)>>s; }
    PlotStream & operator >= ( string*& s)     { return readc(6)>>s; }
    template<class T>    
    PlotStream & operator >= ( KN<T>& b)   { return readc(10), readc(sizeof(T)), operator>>(b);}
    PlotStream & readc(int cc) { int c; read(c); assert(c==cc); return *this;}
    
    void SkipData() { int c; read(c); 
	bool b;
	int i;
	long l,n;
	long long ll;
	string s;
	double d;
	char buf[100];
	switch (c) {
	    case 1: read(b);break;
	    case 2: read(l);break;
	    case 3: read(ll);break;
	    case 4: read(i);break;
	    case 5: read(d);break;
	    case 6: read(s);break;
            case 10: read(l); assert(l>0 && l <100);
		read(n);
		for(int i=0;i<n;++n) read(buf,l);
	    default:
		break;
	}
	
    }    
    
};

#endif //PLOT_STREAM_HPP
