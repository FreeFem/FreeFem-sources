/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARHBMat.h
   Matrix template that generates a matrix in CSC format
   from a Harwell-Boing matrix file.

   ARPACK authors:
      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#ifndef ARHBMAT_H
#define ARHBMAT_H

#include <stddef.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include "arch.h"
#include "arerror.h"


template<class INT, class TYPE>
class ARhbMatrix {

 private:

  char* datafile;        // Filename.
  char  title[73];       // Title.
  char  name[9];         // Name.
  char  type[4];         // Matrix type.
  int   m;               // Number of rows.
  int   n;               // Number of columns.
  int   nnz;             // Number of nonzero variables.
  INT*  irow;            // Row indices.
  INT*  pcol;            // Column pointers.
  TYPE* val;             // Numerical values of matrix entries.

  void ConvertDouble(char* num);

  bool ReadEntry(ifstream& file, int nval, int fval, int& j, double& val);

  bool ReadEntry(ifstream& file, int nval, int fval, int& j, float& val);

  bool ReadEntry(ifstream& file, int nval, int fval,
                 int& j, arcomplex<double>& val);

  bool ReadEntry(ifstream& file, int nval, int fval,
                 int& j, arcomplex<float>& val);

  void ReadFormat(ifstream& file, int& n, int& fmt);

 public:

  bool IsDefined() { return (m!=0); }
  
  bool IsReal() { return (type[0]=='R'); }

  bool IsComplex() { return (type[0]=='C'); }

  bool IsSymmetric() { return (type[1]=='S'); }

  bool IsUnsymmetric() { return (type[1]=='U'); }

  bool IsHermitian() { return (type[1]=='H'); }

  bool IsSkewSymmetric() { return (type[1]=='Z'); }

  char* Filename() { return datafile; }

  char* Title() { return title; }

  char* Name() { return name; }

  char* Type() { return type; }

  int NRows() { return m; }

  int NCols() { return n; }

  int NonZeros() { return nnz; }

  INT* RowInd() { return irow; }

  INT* ColPtr() { return pcol; }

  TYPE* Entries() { return val; }

  void Define(char* filename);
  // Function that reads the matrix file. 

  ARhbMatrix();
  // Short constructor.

  ARhbMatrix(char* filename) { Define(filename); }
  // Long constructor.

  ~ARhbMatrix();
  // Destructor.

}; // Class ARhbMatrix.


// ------------------------------------------------------------------------ //
// ARhbMatrix member functions definition.                                  //
// ------------------------------------------------------------------------ //


template<class INT, class TYPE>
inline void ARhbMatrix<INT, TYPE>::ConvertDouble(char* num) 
{

  char* pd;

  pd = strchr((char*)num,'D');
  if (pd) *pd = 'E';
  pd = strchr((char*)num,'d');
  if (pd) *pd = 'E';


} // ConvertDouble.


template<class INT, class TYPE>
inline bool ARhbMatrix<INT, TYPE>::
ReadEntry(ifstream& file, int nval, int fval, int& j, double& val)
{

  char num[81];
  char c;

  if (file.get((char*)num,fval,'\n')) {
    ConvertDouble((char*)num);
    val = atof((char*)num);
    if (!((++j)%nval)) do file.get(c); while (c!='\n'); 
    return true;
  }
  else {
    return false;
  }

} // ReadEntry (double).


template<class INT, class TYPE>
inline bool ARhbMatrix<INT, TYPE>::
ReadEntry(ifstream& file, int nval, int fval, int& j, float& val)
{

  double dval;
  bool   ret;
  
  ret = ReadEntry(file, nval, fval, j, dval);
  val = (float)dval;
  return ret;

} // ReadEntry (float).


template<class INT, class TYPE>
inline bool ARhbMatrix<INT, TYPE>::
ReadEntry(ifstream& file, int nval, int fval,
          int& j, arcomplex<double>& val)
{

  char num[81], img[81];
  char c;

  if (file.get((char*)num,fval,'\n')) {
    ConvertDouble((char*)num);
    if (!((++j)%nval)) do file.get(c); while (c!='\n'); 
    if (file.get((char*)img,fval,'\n')) {
      ConvertDouble((char*)img);
      if (!((++j)%nval)) do file.get(c); while (c!='\n'); 
      val = arcomplex<double>(atof((char*)num), atof((char*)img));
      return true;
    }
    else {
      return false;
    }
  }
  else {
    return false;
  }

} // ReadEntry (arcomplex<double>).


template<class INT, class TYPE>
inline bool ARhbMatrix<INT, TYPE>::
ReadEntry(ifstream& file, int nval, int fval,
          int& j, arcomplex<float>& val)
{

  // I hope one day c++ will have a standard complex
  // class, so functions like this can be suppressed.

  char num[81], img[81];
  char c;

  if (file.get((char*)num,fval,'\n')) {
    ConvertDouble((char*)num);
    if (!((++j)%nval)) do file.get(c); while (c!='\n'); 
    if (file.get((char*)img,fval,'\n')) {
      ConvertDouble((char*)img);
      if (!((++j)%nval)) do file.get(c); while (c!='\n'); 
      val = arcomplex<float>(atof((char*)num), atof((char*)img));
      return true;
    }
    else {
      return false;
    }
  }
  else {
    return false;
  }

} // ReadEntry (arcomplex<float>).


template<class INT, class TYPE>
void ARhbMatrix<INT, TYPE>::ReadFormat(ifstream& file, int& n, int& fmt)
{

  char c;

  do file.get(c); while ((c != '(') && (c!='\n'));
  file >> n;
  file.get(c);
  while ((c!='I') && (c!='i') && (c!='E') && (c!='e') &&
         (c!='D') && (c!='d') && (c!='\n')) {
    do file.get(c); while ((c != ',') && (c!='\n'));  
    file >> n;
    file.get(c);
  }
  if ((c==')')||(c=='\n')) { // Reading error!
    fmt = 0;
  }
  else {
    file >> fmt;
  }

} // ReadFormat.


template<class INT, class TYPE>
void ARhbMatrix<INT, TYPE>::Define(char* filename)
{

  // Declaring variables.

  int  i, j;
  int  lintot, linptr, linind, linval, linrhs; 
  int  npcol, fpcol, nirow, firow, nval, fval;
  char c;
  char num[81];
  TYPE value;

  // Opening file.

  datafile = filename;
  ifstream file(datafile);
  
  if (!file) {
    throw ArpackError(ArpackError::CANNOT_OPEN_FILE, "ARhbMatrix");
  }

  // Reading the first line.

  file.get((char*)title,73,'\n');
  file.get((char*)name,9,'\n');
  do file.get(c); while (c!='\n'); 

  // Reading the second line.

  file >> lintot >> linptr >> linind >> linval >> linrhs;
  do file.get(c); while (c!='\n'); 

  if ((linptr < 1) || (linind < 1)) { 
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARhbMatrix");
  }

  // Reading the third line.

  file.get((char*)type,4,'\n');
  file >> m >> n >> nnz;
  do file.get(c); while (c!='\n'); 

  if (((type[0] != 'R') && (type[0] != 'C')) || (type[2] != 'A')) {
    throw ArpackError(ArpackError::WRONG_MATRIX_TYPE, "ARhbMatrix");
  }
  else if ((m < 1) || (n < 1) || (nnz < 1)) {
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARhbMatrix");
  }

  // Reading the fourth line.

  ReadFormat(file, npcol, fpcol);
  ReadFormat(file, nirow, firow);
  ReadFormat(file, nval, fval);
  do file.get(c); while (c!='\n'); 
  if ((fpcol<1) || (firow<1) || (fval<1)) {
    throw ArpackError(ArpackError::WRONG_DATA_TYPE, "ARhbMatrix");
  }

  // Skipping the fifth line.

  if (linrhs) {
    do file.get(c); while (c!='\n'); 
    ArpackError(ArpackError::RHS_IGNORED, "ARhbMatrix");
  }

  // Reading column pointers.

  pcol = new INT[n+1];
  fpcol++;
  i = 0;
  while ((i <= n) && (file.get((char*)num,fpcol,'\n'))) {
    pcol[i++] = atoi((char*)num)-1;
    if (!(i%npcol)) do file.get(c); while (c!='\n'); 
  }
  if (i%npcol) do file.get(c); while (c!='\n'); 

  if (i <= n) {
    throw ArpackError(ArpackError::UNEXPECTED_EOF, "ARhbMatrix");
  }

  // Reading row indices.

  irow = new INT[nnz];
  firow++;
  i = 0;
  while ((i < nnz) && (file.get((char*)num,firow,'\n'))) {
    irow[i++] = atoi((char*)num)-1;
    if (!(i%nirow)) do file.get(c); while (c!='\n'); 
  }
  if (i%nirow) do file.get(c); while (c!='\n'); 
  
  if (i < nnz) { 
    throw ArpackError(ArpackError::UNEXPECTED_EOF, "ARhbMatrix");
  }

  // Reading matrix elements.

  fval++;
  val = new TYPE[nnz];
  i = 0;
  j = 0;
  while ((i < nnz) && (ReadEntry(file, nval, fval, j, value))) {
    val[i++] = value;
  }  
  if (j%nval) do file.get(c); while (c!='\n'); 

  if (i < nnz) {
    throw ArpackError(ArpackError::UNEXPECTED_EOF, "ARhbMatrix");
  }

  // Closing file and reporting success.

  file.close();

} // Define.


template<class INT, class TYPE>
ARhbMatrix<INT, TYPE>::ARhbMatrix()
{

  m = n = nnz = 0;
  title[0]= '\0';
  name[0] = '\0';
  type[0] = '\0';
  pcol    = NULL;
  irow    = NULL;
  val     = NULL;

} // Short constructor.


template<class INT, class TYPE>
ARhbMatrix<INT, TYPE>::~ARhbMatrix()
{

  if (irow != NULL) delete[] irow;
  if (pcol != NULL) delete[] pcol;
  if (val  != NULL) delete[] val;

} // Destructor.


#endif ARHBMAT_H

