/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARGNSym.h.
   Arpack++ class ARNonSymGenEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARGNSYM_H
#define ARGNSYM_H

#include <stddef.h>
#include "arch.h"
#include "blas1c.h"
#include "lapackc.h"
#include "arsnsym.h"
#include "argeig.h"
#include "arrgnsym.h"

template<class FLOAT, class FOP, class FB>
class ARNonSymGenEig:
  virtual public ARGenEig<FLOAT, FLOAT, FOP, FB>,
  virtual public ARNonSymStdEig<FLOAT, FOP>,
  virtual public ARrcNonSymGenEig<FLOAT>  {

 protected:

 // a) Protected variables:

  FB     *objA;      // Object that has MultAx as a member function.
  TypeBx MultAx;     // Function that evaluates the product A*x.


 // b) Protected functions:

  void RecoverEigenvalues();
  // Uses Rayleigh quotient to recover eigenvalues of the original
  // problem when shift is complex.

  virtual void Copy(const ARNonSymGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  virtual void SetShiftInvertMode(FLOAT sigmaRp, FOP* objOPp,
                                  void (FOP::* MultOPxp)(FLOAT[], FLOAT[]));
  // Turns the problem to real shift-and-invert mode with sigmaRp as shift.

  virtual void SetComplexShiftMode(char partp, FLOAT sigmaRp, 
                                   FLOAT sigmaIp, FOP* objOPp, 
                                   void (FOP::* MultOPxp)(FLOAT[], FLOAT[]), 
                                   FB* objAp,
                                   void (FB::* MultAxp)(FLOAT[], FLOAT[]));
  // Turns the problem to complex shift-and-invert mode with shift
  // defined by sigmaRp and sigmaIp. MultAx is used to obtain eigenvalues.


 // c.2) Functions that perform all calculations in one step.

  virtual int FindEigenvalues();
  // Determines nev approximated eigenvalues of the given eigen-problem.

  virtual int FindEigenvectors(bool schurp = false);
  // Determines nev approximated eigenvectors of the given eigen-problem
  // Optionally also determines nev Schur vectors that span the desired
  // invariant subspace.

  virtual int FindSchurVectors();
  // Determines nev Schur vectors that span the desired invariant subspace.
  // Redefined in ARSymEig.


 // c.3) Constructors and destructor.

  ARNonSymGenEig() { part = 'R'; }
  // Short constructor (Does nothing but calling base classes constructors).

  ARNonSymGenEig(int np, int nevp, FOP* objOPp,
                 void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
                 FB* objBp, void (FB::* MultBxp)(FLOAT[], FLOAT[]),
                 char* whichp = "LM", int ncvp = 0, FLOAT tolp = 0.0,
                 int maxitp = 0, FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARNonSymGenEig(int np, int nevp, FOP* objOPp,
                 void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
                 FB* objBp, void (FB::* MultBxp)(FLOAT[], FLOAT[]),
                 FLOAT sigmap, char* whichp = "LM", int ncvp = 0,
                 FLOAT tolp = 0.0, int maxitp = 0, FLOAT* residp = NULL,
                 bool ishiftp = true);
  // Long constructor (real shift and invert mode).

  ARNonSymGenEig(int np, int nevp, FOP* objOPp,
                 void (FOP::* MultOPxp)(FLOAT[], FLOAT[]), FB* objAp,
                 void (FB::* MultAxp)(FLOAT[], FLOAT[]), FB* objBp,
                 void (FB::* MultBxp)(FLOAT[], FLOAT[]), char partp,
                 FLOAT sigmaRp, FLOAT sigmaIp, char* whichp = "LM",
                 int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                 FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (complex shift and invert mode).

  ARNonSymGenEig(const ARNonSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARNonSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARNonSymGenEig& operator=(const ARNonSymGenEig& other);
  // Assignment operator.

}; // class ARNonSymGenEig.


// ------------------------------------------------------------------------ //
// ARNonSymGenEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class FLOAT, class FOP, class FB>
inline void ARNonSymGenEig<FLOAT, FOP, FB>::
Copy(const ARNonSymGenEig<FLOAT, FOP, FB>& other)
{

  ARGenEig<FLOAT, FLOAT, FOP, FB>::Copy(other);
  objA   = other.objA;
  MultAx = other.MultAx;
  part   = other.part;

} // Copy.


template<class FLOAT, class FOP, class FB>
void ARNonSymGenEig<FLOAT, FOP, FB>::RecoverEigenvalues()
{

  int    j, ColJ, ColJp1;
  FLOAT  numr, numi, denr, deni;
  FLOAT* Ax;

  Ax = new FLOAT[n];

  for (j=0; j<nconv; j++) {

    ColJ   = j*n;
    ColJp1 = ColJ+n;

    if (EigValI[j] == (FLOAT)0.0) {

      // Eigenvalue is real. Computing EigVal = x'(Ax)/x'(Mx).

      (objB->*MultAx)(&EigVec[ColJ], Ax);
      numr = dot(n, &EigVec[ColJ], 1, Ax, 1);
      (objB->*MultBx)(&EigVec[ColJ], Ax);
      denr = dot(n, &EigVec[ColJ], 1, Ax, 1);
      EigValR[j] =  numr / denr;

    }
    else {

      // Eigenvalue is complex.

      // Computing x'(Ax).

      (objB->*MultAx)(&EigVec[ColJ], Ax);
      numr = dot(n, &EigVec[ColJ], 1, Ax, 1);
      numi = dot(n, &EigVec[ColJp1], 1, Ax, 1);
      (objB->*MultAx)(&EigVec[ColJp1], Ax);
      numr = numr + dot(n, &EigVec[ColJp1], 1, Ax, 1);
      numi = -numi + dot(n, &EigVec[ColJ], 1, Ax, 1);

      // Computing x'(Mx).

      (objB->*MultBx)(&EigVec[ColJ], Ax);
      denr = dot(n, &EigVec[ColJ], 1, Ax, 1);
      deni = dot(n, &EigVec[ColJp1], 1, Ax, 1);
      (objB->*MultBx)(&EigVec[ColJp1], Ax);
      denr = denr + dot(n, &EigVec[ColJp1], 1, Ax, 1);
      deni = -deni + dot(n, &EigVec[ColJ], 1, Ax, 1);

      // Computing the first eigenvalue of the conjugate pair.

      EigValR[j] = (numr*denr+numi*deni) / lapy2(denr, deni);
      EigValI[j] = (numi*denr-numr*deni) / lapy2(denr, deni);

      // Getting the second eigenvalue of the conjugate pair by taking
      // the conjugate of the first.

      EigValR[j+1] = EigValR[j];
      EigValI[j+1] = -EigValI[j];
      j++;

    }

  }

  delete[] Ax;

} // RecoverEigenvalues.


template<class FLOAT, class FOP, class FB>
inline void ARNonSymGenEig<FLOAT, FOP, FB>::
SetShiftInvertMode(FLOAT sigmaRp, FOP* objOPp,
                   void (FOP::* MultOPxp)(FLOAT[], FLOAT[]))
{

  part    = 'R';
  objOP   = objOPp;
  MultOPx = MultOPxp;
  ChangeShift(sigmaRp);

} // SetShiftInvertMode.


template<class FLOAT, class FOP, class FB>
inline void ARNonSymGenEig<FLOAT, FOP, FB>::
SetComplexShiftMode(char partp, FLOAT sigmaRp, FLOAT sigmaIp, 
                    FOP* objOPp, void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
                    FB* objAp, void (FB::* MultAxp)(FLOAT[], FLOAT[]))
{

  objOP   = objOPp;
  MultOPx = MultOPxp;
  objA    = objAp;
  MultAx  = MultAxp;
  part    = CheckPart(partp);
  ChangeShift(sigmaRp, sigmaIp);

} // SetComplexShiftMode.


template<class FLOAT, class FOP, class FB>
inline int ARNonSymGenEig<FLOAT, FOP, FB>::FindEigenvalues()
{

  nconv = ARStdEig<FLOAT, FLOAT, FOP>::FindEigenvalues();
  if (sigmaI != 0.0) RecoverEigenvalues();
  return nconv;

} // FindEigenvalues.


template<class FLOAT, class FOP, class FB>
inline int ARNonSymGenEig<FLOAT, FOP, FB>::FindEigenvectors(bool schurp)
{

  nconv = ARStdEig<FLOAT, FLOAT, FOP>::FindEigenvectors(schurp);
  if (sigmaI != 0.0) RecoverEigenvalues();
  return nconv;

} // FindEigenvectors.


template<class FLOAT, class FOP, class FB>
int ARNonSymGenEig<FLOAT, FOP, FB>::FindSchurVectors()
{

  nconv = ARStdEig<FLOAT, FLOAT, FOP>::FindSchurVectors();
  if (sigmaI != 0.0) RecoverEigenvalues();
  return nconv;

} // FindSchurVectors.


template<class FLOAT, class FOP, class FB>
inline ARNonSymGenEig<FLOAT, FOP, FB>::
ARNonSymGenEig(int np, int nevp, FOP* objOPp,
               void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
               FB* objBp, void (FB::* MultBxp)(FLOAT[], FLOAT[]),
               char* whichp, int ncvp, FLOAT tolp, int maxitp,
               FLOAT* residp, bool ishiftp)

{

  part = 'R';                // Considering mode = 3 in ChangeShift.
  NoShift();
  DefineParameters(np, nevp, objOPp, MultOPxp, objBp, MultBxp,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT, class FOP, class FB>
inline ARNonSymGenEig<FLOAT, FOP, FB>::
ARNonSymGenEig(int np, int nevp, FOP* objOPp,
               void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
               FB* objBp, void (FB::* MultBxp)(FLOAT[], FLOAT[]),
               FLOAT sigmap, char* whichp, int ncvp,
               FLOAT tolp, int maxitp, FLOAT* residp, bool ishiftp)

{

  SetShiftInvertMode(sigmap, objOPp, MultOPxp);
  DefineParameters(np, nevp, objOPp, MultOPxp, objBp, MultBxp,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);


} // Long constructor (real shift and invert mode).


template<class FLOAT, class FOP, class FB>
inline ARNonSymGenEig<FLOAT, FOP, FB>::
ARNonSymGenEig(int np, int nevp, FOP* objOPp,
               void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
               FB* objAp, void (FB::* MultAxp)(FLOAT[], FLOAT[]),
               FB* objBp, void (FB::* MultBxp)(FLOAT[], FLOAT[]),
               char partp, FLOAT sigmaRp, FLOAT sigmaIp,
               char* whichp, int ncvp, FLOAT tolp, int maxitp,
               FLOAT* residp, bool ishiftp)

{

  SetComplexShiftMode(partp, sigmaRp, sigmaIp, objOPp,
                      MultOPxp, objAp, MultAxp);
  DefineParameters(np, nevp, objOPp, MultOPxp, objBp, MultBxp,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class FLOAT, class FOP, class FB>
ARNonSymGenEig<FLOAT, FOP, FB>& ARNonSymGenEig<FLOAT, FOP, FB>::
operator=(const ARNonSymGenEig<FLOAT, FOP, FB>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARGNSYM_H

