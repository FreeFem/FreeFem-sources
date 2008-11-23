/*
 *  PlotStream.hpp
 *  ff
 *
 *  Created by FrÃ©dÃ©ric Hecht on 27/10/08.
 *  Copyright 2008 UPMC. All rights reserved.
 *
 */
#include <cstdio>

using  Fem2D::Mesh;
class PlotStream 
{
public:
  
  FILE * TheStream; 
  PlotStream(FILE *Stream) :TheStream(Stream) {}
  operator bool() const { return TheStream;}
  // datatype mush be < 0 to have no collistion with arg number. 
  enum datatype { dt_meshes=-1,dt_plots=-2,dt_endplot=-3,dt_endarg=99999,dt_newplot=-5  };
  
  void SendNewPlot() {  write((long )dt_newplot); }
  void SendEndArgPlot() {write((long )dt_endarg); }
  void SendEndPlot() { write((long )dt_endplot);fflush(TheStream); }
  void SendPlots() { write((long )dt_plots); }
  void SendMeshes() { write((long )dt_meshes);}
  void write(const void *data,size_t l) {fwrite(data,1,l,TheStream);}
  PlotStream& write(const bool& b) {write(reinterpret_cast<const void *> (&b),sizeof(bool));return *this;}
  PlotStream& write(const long& b) {write(reinterpret_cast<const void *> (&b),sizeof(long));return *this;}
  PlotStream& write(const double& b) {write(reinterpret_cast<const void *> (&b),sizeof(double));return *this;}
  PlotStream& write(const string& b) {  
    size_t l=b.size();
    // cout << " l : " << b.size() <<endl; 
    write(reinterpret_cast<const void *> (&l),sizeof(size_t));
    write(b.data(),l);
    return *this;
    }
  
  
    
  PlotStream & operator << (const bool& b) 
  { return write(b); }
  PlotStream & operator << (const long& b) 
  { return write(b); }        
  PlotStream & operator << (const double& b) 
  { return write(b); }
  PlotStream & operator << (const string& s) 
  { return write(s); }
  PlotStream & operator << (const string* s) 
  { ffassert(s); return write(*s); }
  PlotStream & operator << (const Mesh& Th) {
    Serialize s=Th.serialize();
    long  n=s.size();
    write( n );
    write(s,s.size());
    return *this;}
  PlotStream & operator << (const KN_<double>& b)
  {
    long n=b.N();
    write(n);
    if(b.end() != n)
      {
	for (int i=0;i<n;++i)
	  write(b[i]);
      }
    else 
      write(reinterpret_cast<const char *> ((double *) b),n*sizeof(double));  
    return *this;
  }
  
  void read( void *data,size_t l) {
	char * p= (char*)data;
	for(int i=0;i<l;++i)	
	  *p++ = (char) getc(TheStream);
//	fread(  p++,1,1,TheStream);
//       read(data,l);
}
  bool good() const {return ferror(TheStream)==0;}
  void GetNewPlot() { get(dt_newplot) ;}
  void GetEndArgPlot() {get(dt_endarg); }
  void GetEndPlot() {get(dt_endplot); }
  void GetPlots() { get(dt_plots); }
  void GetMeshes() { get(dt_meshes);}
  void get(datatype t) { long tt; read(tt);
    if( tt !=(long) t) cout << " Error Check :  get " << tt << " == wiat for  "<< t << endl; 
    ffassert(tt==(long) t);}
  PlotStream& read( bool& b) {read(reinterpret_cast< void *> (&b),sizeof(bool));return *this;}
  PlotStream& read( long& b) {read(reinterpret_cast< void *> (&b),sizeof(long));return *this;}
  PlotStream& read( double& b) {read(reinterpret_cast< void *> (&b),sizeof(double));return *this;}
  PlotStream& read( string& b) {  	
    size_t l=b.size();
    read(reinterpret_cast< void *> (&l),sizeof(size_t));
    b.resize(l);
    //    cout << "read str :len str=" << l << endl;
    read(& (b[0]),l);
    return *this;
  }
  
  
  
  PlotStream & operator >> ( bool& b) 
  { return read(b); }
  PlotStream & operator >> ( long& b) 
  { return read(b); }        
  PlotStream & operator >> ( double& b) 
  { return read(b); }
  PlotStream & operator >> ( string& s) 
  { return read(s); }
  PlotStream & operator >> ( string *& s) 
  { if(!s) s= new string; return read(*s);
    // cout << " fread string " << s <<endl;
  }
  
  PlotStream & operator >> ( Mesh *& Th) {
    long n;
    read(n);
    Serialize s(n,Mesh::magicmesh);
    read(s,n );
    Th= new Mesh(s);
    return *this;}
  
  PlotStream & operator >> ( KN<double>& b)
  {
    long n;
    read(n);
    if( ! b.N() ) b.init(n);
    ffassert( b.N()==n); 
    if(b.end() != n)
      {
	for (int i=0;i<n;++i)
	  read(b[i]);
	  }
    else 
      read(reinterpret_cast< char *> ((double *) b),n*sizeof(double));
    //     cout << "PlotStream read : " << n << " " << b.min() << " " << b.max() << endl;
    return *this;
  }
  
  
};

