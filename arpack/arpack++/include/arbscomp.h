/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARBSComp.h.
   Arpack++ class ARluCompStdEig definition
   (band matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARBSCOMP_H
#define ARBSCOMP_H

#include <stddef.h>
#include "arch.h"
#include "arscomp.h"
#include "arbnsmat.h"
#include "arrseig.h"


template<class FLOAT>
class ARluCompStdEig:
  public virtual ARCompStdEig<FLOAT, ARbdNonSymMatrix<arcomplex<FLOAT> > > {

 public:

 // a) Public functions:

 // a.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(arcomplex<FLOAT> sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(arcomplex<FLOAT> sigmap);

 // a.2) Constructors and destructor.

  ARluCompStdEig() { }
  // Short constructor.

  ARluCompStdEig(int nevp, ARbdNonSymMatrix<arcomplex<FLOAT> >& A, 
                 char* whichp = "LM", int ncvp = 0,
                 FLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<FLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluCompStdEig(int nevp, ARbdNonSymMatrix<arcomplex<FLOAT> >& A,
                 arcomplex<FLOAT> sigma, char* whichp = "LM",
                 int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<FLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARluCompStdEig(const ARluCompStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluCompStdEig() { }
  // Destructor.


 // b) Operators.

  ARluCompStdEig& operator=(const ARluCompStdEig& other);
  // Assignment operator.

}; // class ARluCompStdEig.


// ------------------------------------------------------------------------ //
// ARluCompStdEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class FLOAT>
inline void ARluCompStdEig<FLOAT>::
ChangeShift(arcomplex<FLOAT> sigmaRp)
{

  objOP->FactorAsI(sigmaRp);
  ARrcStdEig<FLOAT, arcomplex<FLOAT> >::ChangeShift(sigmaRp);

} // ChangeShift.


template<class FLOAT>
inline void ARluCompStdEig<FLOAT>::SetRegularMode()
{

  ARStdEig<FLOAT, arcomplex<FLOAT>, ARbdNonSymMatrix<arcomplex<FLOAT> > >::
    SetRegularMode(objOP, &ARbdNonSymMatrix<arcomplex<FLOAT> >::MultMv);

} // SetRegularMode.


template<class FLOAT>
inline void ARluCompStdEig<FLOAT>::
SetShiftInvertMode(arcomplex<FLOAT> sigmap)
{

  ARStdEig<FLOAT, arcomplex<FLOAT>, ARbdNonSymMatrix<arcomplex<FLOAT> > >::
    SetShiftInvertMode(sigmap, objOP,
                       &ARbdNonSymMatrix<arcomplex<FLOAT> >::MultInvv);

} // SetShiftInvertMode.


template<class FLOAT>
inline ARluCompStdEig<FLOAT>::
ARluCompStdEig(int nevp, ARbdNonSymMatrix<arcomplex<FLOAT> >& A,
               char* whichp, int ncvp, FLOAT tolp,
               int maxitp, arcomplex<FLOAT>* residp, bool ishiftp)

{

  NoShift();
  DefineParameters(A.ncols(), nevp, &A,
                   ARbdNonSymMatrix<arcomplex<FLOAT> >::MultMv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARluCompStdEig<FLOAT>::
ARluCompStdEig(int nevp, ARbdNonSymMatrix<arcomplex<FLOAT> >& A,
               arcomplex<FLOAT> sigmap, char* whichp, int ncvp,
               FLOAT tolp, int maxitp, arcomplex<FLOAT>* residp,
               bool ishiftp)

{

  DefineParameters(A.ncols(), nevp, &A, 
                   &ARbdNonSymMatrix<arcomplex<FLOAT> >::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class FLOAT>
ARluCompStdEig<FLOAT>& ARluCompStdEig<FLOAT>::
operator=(const ARluCompStdEig<FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBSCOMP_H
