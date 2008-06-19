// -*- Mode : c++ -*-
//
// SUMMARY  :      
// USAGE    :        
// ORG      : 
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

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

#include "meshtype.h"
#include "Meshio.h"
namespace bamg {

const char *OFortranUnFormattedFile::unkown("unkown file name");
void (*MeshIstreamErrorHandler)(ios & ) =0;
///////////////////////////////////////////////////////////
 void WriteStr(ostream & out,char * str)
{ 
  int i=0;
  char c;
  out << '"' ;
  while (i++,c=*str++) {
   if (i==70) 
      cout << " //\n",i=0;
   if (c == '"') out << "\"\"" ;
   else out << c ;}
 out << '"' ;
}
///////////////////////////////////////////////////////////
double * ReadbbFile(const char * file,long & nbsol,long & lsol,const int typesol,const int dim)
{

  MeshIstream frbb(file);
  int dimlu,typesollu;
  long i;
 double *sol =0;
  frbb >>    dimlu >> nbsol >> lsol >> typesollu ;
  if(typesol != typesollu )
    { 
      cerr << " incorrect type of solution (read) " << typesollu << " != (wanted) " <<  typesol << endl;
      cerr << "  or       dim  of solution (read) " << dimlu  << " != (wanted) " <<  dim << endl;
      nbsol=0;
      lsol=0;
    }
  else
    {
      sol = new double[lsol*nbsol];
      double *s=sol;
      for(  i=0;i<lsol;i++)
	for (long j=0;j<nbsol;j++)
	  frbb >>  *s++ ;
      //      frbb.close();
    }

  
  return sol;
}

///////////////////////////////////////////////////////////
double * ReadBBFile(const char * file,long & nbsol,long & lsol,int *& typesols,
  const int dim,const int typesol)
{

  MeshIstream frbb(file);
  int dimlu,typesollu;
  typesols = 0;
  long i,n=0;
  double *sol =0;
  int sizeoftype[] = { 1, dim ,dim * (dim+1) / 2, dim * dim } ;
  // char * charoftype[] = {  "Scalar" , "dim Vector" ,"dim x dim Sym-Matric","dim x dim Full-Matric" } ;
  frbb >>    dimlu >> nbsol ;
  typesols = new int [nbsol] ;
  for (i=0;i<nbsol;i++)
    {
     frbb >> typesols[i];
     typesols[i]--;
     if (typesols[i]<0 || typesols[i]>= 4)
       {
	 cerr << " Error reading BBFile the type solution " << i+1 << " is " << typesols[i]+1 
	      << " is not in [1..4] " << endl;
	   frbb.ShowIoErr(998);
	 nbsol=0;
	 lsol=0;
	 delete [] typesols;
	 typesols=0;
	 return 0; 
       }
     n += sizeoftype[ typesols[i]];
    }
  frbb >> lsol >> typesollu ;
  if(typesol != typesollu )
    { 
      cerr << " incorrect type of solution (read) " << typesollu << " != (wanted) " <<  typesol << endl;
      cerr << "  or       dim  of solution (read) " << dimlu  << " != (wanted) " <<  dim << endl;
      frbb.ShowIoErr(999);
      nbsol=0;
      lsol=0;
    }
  else
    {
      if (verbosity> 5)
	cout << "     read BB file " << file << "  with nbsol " << nbsol << " total nb of field = " << n << endl; 
      if (verbosity>9) 
	for (i=0;i<nbsol;i++)
	  cout << "      the type of solution " << i+1 << " is " <<  sizeoftype[ typesols[i]]  
	       << " and  the number of sub-field are " << sizeoftype[ typesols[i]] << endl; 
      sol = new double[lsol*n];
      double *s=sol;
      for(  i=0;i<lsol;i++)
	for (long j=0;j<n;j++)
	  frbb >>  *s++ ;
      //      frbb.close();
    }

  
  return sol;
}


void MeshIstream::ShowIoErr(int s) {
  LineError = 1;
  if (CurrentFile)
    cerr << " In  file " <<  CurrentFile ;
  cerr << ", Erreur Lecture "  << s << 
    ", good="  << (s & ios::goodbit) <<
    ", bad=" << (s & ios::badbit) <<
    ", fail=" << (s & ios::failbit) <<
    ", eof=" << (s & ios::eofbit) << " Line " << LineNumber <<
    endl;
  if(!in.eof()) 
    {
      in.clear(ios::goodbit);
      int i=0;
      char c;
      cerr << " the next character :" ;
      while (i++<80 && in.get(c) && (c != '\n' || i < 30 ))
	cerr.put(c);
      cerr << endl; 
    }
  if (MeshIstreamErrorHandler)
    MeshIstreamErrorHandler(in);
  else 
   { 
    in.clear(ios::failbit);
   }
}
int  MeshIstream::IsString(const char* s)
{ 
  int n=0;
  char c;
  const char * ss = s;
  while (in.get(c) && c ==' ') n++; // eat whitespace 
  if (in.good())
    while ( *ss && c == *ss && in.get(c) )
      ss++; 

  if (*ss) { // no 
    if (verbosity>9) 
      cout << "IsString: not " << s << " " << n  << "  putback " << ss-s << " :" << c  ;
    if (in.good())
      in.putback(c),ss--;
    while(ss-s>0)
      {
	cout << *ss ;
	in.putback(*ss--);
      };
    if (verbosity>9) 
      cout << ":"<< endl;
    return 0;}
  else return 1;
}

 char * MeshIstream::ReadStr()
{ 
  static char buf[1025];
  int  instr=0,k=0;
  char * b=buf, bb=0, *bend=buf+1023;
  
  for (b=buf;b<bend &&  in.get(*b);instr && b++)
    {
      int kold = k;
      if(isspace(*b))
	{
	  if (*b == '\n') LineNumber++;
	  if (instr && !bb) break;
	}
      else if (*b == bb) 
	{
	  if(++k%2) b--;
	}
      else if( !instr) 
	{
	  instr =1;
	  if  ((*b == '"') ||  (*b == '\'' ) )
	    bb=*b--,k=0;
	}
      else
	instr =1;
      
      if( k%2 && kold==k ) // even and fin
	break;
    }
  *b=0;
  in.clear(ios::goodbit);
  //  cout << " l = " << b-buf << " " ;
  char * r = new char [b-buf+1];
  strcpy(r,buf);
  return r ;  
}
///////////////////////////////////////////////////////////
OFortranUnFormattedFile::~OFortranUnFormattedFile()
{
  // save the last record 
  if (l==0) { // l rec no set --
    l=where();
  if (verbosity>9)
    cout << " size of last record  = " 
	 << l << " n = " << nb_rec << " i= " << i << endl;
   if (!f->good()) Error(3);
    f->seekp(i-sizeof(long));
   if (!f->good()) Error(3);
    f->write((char*)&l,sizeof(l)); 
    if (!f->good()) Error(3);
   n=j;
  }

  f->seekp(n);
  if (!f->good()) Error(3);
   f->write((char*) &l,sizeof(l));
  if (!f->good()) Error(3);
   
  if (f && to_close) {
    if (verbosity>9)
      cout << "delete OFortranUnFormattedFile " << file_name 
	   << " @end  = " << n <<  endl;
    delete f;} f=0;
}

IFortranUnFormattedFile::~IFortranUnFormattedFile()
{
  if (f && to_close) {
    if (verbosity>9)
      cout << " delete IFortranUnFormattedFile" << file_name 
	   << "  @end = " << n <<  endl;
    delete f;} f=0;
}

long IFortranUnFormattedFile::Record()
{
  nb_rec++;
  n += sizeof(long);
  f->seekg(n);
  f->read((char*)&l,sizeof(l));
  if (!f->good()) Error(3);
  i=j= n+sizeof(l); 
  n +=   l + sizeof(l) ; // end 
  if (verbosity>9)
  cout << " Read  rec end =" << n << " l=  " << l 
       << " current= " << j << " begin= " << i << endl;
  return l;
}

long OFortranUnFormattedFile::Record(long ll) //
{
  if (j==4 && l==0) 
    { 
      l=ll;
      f->seekp(0);
      f->write((char*)&ll,sizeof(l));
      if (!f->good()) Error(3);
      return 0;
    }
  if (n>=0) 
    {
      if (l==0) 
	{ // l rec no set --
	  l=where();
	  if (verbosity>9)
	    cout << " set len  of rec  " << nb_rec  << " = " << l 
		 << " @ " << i-sizeof(long) << endl;

	  f->seekp(i-sizeof(long));
	  f->write((char*)&l,sizeof(l)); 
	  n=j;
	}
      f->seekp(n);
      f->write((char*)&l,sizeof(l)); 
      if (!f->good()) Error(3);
    }
  else 
    f->seekp(0);
  n += sizeof(long);
  nb_rec++;
  f->write((char*)&ll,sizeof(l));
  if (!f->good()) Error(3);
  
  l=ll;
  n +=  sizeof(long);
  j = n;
  i = n;
  n +=   l;
  if (verbosity>9)
    cout << " Write rec end =" << n << " l=  " << l 
	 << " current= " << j << " begin= " << i << endl;
  return l;
}
void OFortranUnFormattedFile::Error(int err)
{
  cerr << "Fatal Error Output FortranUnFormattedFile " << endl;
  if (err==0) 
    cerr << "-- Error openning ";
  else if (err==1)
    cerr<< "-- OverFlow write ";
  else if (err==2)
    cerr<< "-- Erreur  write " ;
  else if (err==3)
    cerr<< "-- Erreur  write record info " ;
  else 
    cerr << " unkown err " << err ;

  cerr << " Record number = "    << nb_rec << endl
       << " read position  in file " << j << " < " << n 
       << " =  end on record " << endl;
  cerr << " position in the record = " << where() 
       << " length of record = " << l <<endl;
  cerr 	<< " file = " << file_name<< endl ;
  MeshError(900);
  //  throw(ErrorMesh("exit",1));
}
void IFortranUnFormattedFile::Error(int err)
{
  cerr << "Fatal Error Input FortranUnFormattedFile " << endl;
  if (err==0) 
    cerr << " Error openning ";
  else if (err==1)
    cerr<< " OverFlow read";
  else if (err==2)
    cerr<< " Erreur  read" ;
  else if (err==3)
    cerr<< " Erreur  read record info " ;
  else 
    cerr << " unkown err " << err ;

  cerr << " Record number = "    << nb_rec << endl
       << " read position  in file " << j << " < " << n 
       << " =  end on record " << endl;
  cerr << " position in the record = " << where() 
       << " length of record = " << l <<endl;
  cerr 	<< " file = " << file_name<< endl ;
    MeshError(900);
    //  throw(ErrorMesh("exit",1));
}

}// end of namespace bamg 
