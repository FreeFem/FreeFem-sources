/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARError.h.
   Definition of ArpackError, a class that handles errors
   occurred during Arpack execution.

   There are three ways of handling an error:
   a) Declaring a variable of type ArpackError and calling
      function Set with the correct ErrorCode (see codes below).
   b) Calling the constructor ArpackError(ErrorCode) to define
      a variable.
   c) Calling ArpackError::Set(ErrorCode) directly.

   If an error occurs, a brief description of the error is
   displayed on the "cerr" stream, unless the variable
   ARPACK_SILENT_MODE is defined.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARERROR_H
#define ARERROR_H

#include "arch.h"
#include <iostream.h>
#include <stdlib.h>

class ArpackError {

 public:

  enum ErrorCode {  // Listing all kinds of errors.

    // Innocuous error type.

    NO_ERRORS           =    0,

    // Errors in parameter definitions.

    PARAMETER_ERROR     = -101,
    N_SMALLER_THAN_2    = -102,
    NEV_OUT_OF_BOUNDS   = -103,
    WHICH_UNDEFINED     = -104,
    PART_UNDEFINED      = -105,
    INVMODE_UNDEFINED   = -106,
    RANGE_ERROR         = -107,

    // Errors in Aupp and Eupp functions.

    LAPACK_ERROR        = -201,
    START_RESID_ZERO    = -202,
    NOT_ACCURATE_EIG    = -203,
    REORDERING_ERROR    = -204,
    ARNOLDI_NOT_BUILD   = -205,
    AUPP_ERROR          = -291,
    EUPP_ERROR          = -292,

    // Errors in main functions.

    CANNOT_PREPARE      = -301,
    CANNOT_FIND_BASIS   = -302,
    CANNOT_FIND_VALUES  = -303,
    CANNOT_FIND_VECTORS = -304,
    CANNOT_FIND_SCHUR   = -305,
    SCHUR_UNDEFINED     = -306,

    // Errors due to incorrect function calling sequence.

    CANNOT_GET_VECTOR   = -401,
    CANNOT_GET_PROD     = -402,
    CANNOT_PUT_VECTOR   = -403,
    PREPARE_NOT_OK      = -404,
    BASIS_NOT_OK        = -405,
    VALUES_NOT_OK       = -406,
    VECTORS_NOT_OK      = -407,
    SCHUR_NOT_OK        = -408,
    RESID_NOT_OK        = -409,

    // Errors in classes that perform LU decompositions.

    MATRIX_IS_SINGULAR  = -501,
    DATA_UNDEFINED      = -502,
    INSUFICIENT_MEMORY  = -503,
    NOT_SQUARE_MATRIX   = -504,
    NOT_FACTORED_MATRIX = -505,
    INCOMPATIBLE_SIZES  = -506,
    DIFFERENT_TRIANGLES = -507,
    INCONSISTENT_DATA   = -508,
    CANNOT_READ_FILE    = -509,

    // Errors in matrix files.

    CANNOT_OPEN_FILE    = -551, 
    WRONG_MATRIX_TYPE   = -552,
    WRONG_DATA_TYPE     = -553,
    RHS_IGNORED         = -554,
    UNEXPECTED_EOF      = -555,

    // Other severe errors.

    NOT_IMPLEMENTED     = -901,
    MEMORY_OVERFLOW     = -902,
    GENERIC_SEVERE      = -999,

    // Warnings.

    NCV_OUT_OF_BOUNDS   =  101,
    MAXIT_NON_POSITIVE  =  102,
    MAX_ITERATIONS      =  201,
    NO_SHIFTS_APPLIED   =  202,
    CHANGING_AUTOSHIFT  =  301,
    DISCARDING_FACTORS  =  401,
    GENERIC_WARNING     =  999

  };

 private:

  static ErrorCode code;

  static void Print(const char* where, const char* message);
  // Writes error messages on cerr stream.

 public:

  static void Set(ErrorCode error, char* where="AREigenProblem");
  // Set error code and write error messages.

  static int Status() { return (int) code; }
  // Returns current value of error code.

  ArpackError(ErrorCode error, char* where="AREigenProblem") {
    Set(error,where);
  }
  // Constructor that set error code.

  ArpackError() { code = NO_ERRORS; };
  // Constructor that does nothing.

};

inline void ArpackError::Print(const char* where, const char* message)
{

#ifndef ARPACK_SILENT_MODE
  cerr << "Arpack error in " << where << "." << endl;
  cerr << "-> " << message << "." << endl;
#endif

} // Print

void ArpackError::Set(ErrorCode error, char* where)
{

  code = error;
  switch (code) {
  case NO_ERRORS          :
    return;
  case NOT_IMPLEMENTED    :
    Print(where, "This function was not implemented yet");
    return;
  case MEMORY_OVERFLOW    :
    Print(where, "Memory overflow");
    return;
  case GENERIC_SEVERE     :
    Print(where, "Severe error");
    return;
  case PARAMETER_ERROR    :
    Print(where, "Some parameters were not correctly defined");
    return;
  case N_SMALLER_THAN_2   :
    Print(where, "'n' must be greater than one");
    return;
  case NEV_OUT_OF_BOUNDS  :
    Print(where, "'nev' is out of bounds");
    return;
  case WHICH_UNDEFINED    :
    Print(where, "'which' was not correctly defined");
    return;
  case PART_UNDEFINED     :
    Print(where, "'part' must be one of 'R' or 'I'");
    return;
  case INVMODE_UNDEFINED  :
    Print(where, "'InvertMode' must be one of 'S' or 'B'");
    return;
  case RANGE_ERROR  :
    Print(where, "Range error");
    return;
  case LAPACK_ERROR       :
    Print(where, "Could not perform LAPACK eigenvalue calculation");
    return;
  case START_RESID_ZERO   :
    Print(where, "Starting vector is zero");
    return;
  case NOT_ACCURATE_EIG   :
    Print(where, "Could not find any eigenvalue to sufficient accuracy");
    return;
  case REORDERING_ERROR   :
    Print(where, "Reordering of Schur form was not possible");
    return;
  case ARNOLDI_NOT_BUILD  :
    Print(where, "Could not build an Arnoldi factorization");
    return;
  case AUPP_ERROR         :
    Print(where, "Error in ARPACK Aupd fortran code");
    return;
  case EUPP_ERROR         :
    Print(where, "Error in ARPACK Eupd fortran code");
    return;
  case CANNOT_PREPARE     :
    Print(where, "Could not correctly define internal variables");
    return;
  case CANNOT_FIND_BASIS  :
    Print(where, "Could not find the Arnoldi basis vectors");
    return;
  case CANNOT_FIND_VALUES :
    Print(where, "Could not find any eigenvalue");
    return;
  case CANNOT_FIND_VECTORS:
    Print(where, "Could not find any eigenvector");
    return;
  case CANNOT_FIND_SCHUR  :
    Print(where, "Could not find any Schur vector");
    return;
  case SCHUR_UNDEFINED    :
    Print(where, "FindEigenvectors must be used instead of FindSchurVectors");
    return;
  case CANNOT_GET_VECTOR  :
    Print(where, "Vector is not already available");
    return;
  case CANNOT_GET_PROD    :
    Print(where, "Matrix vector product is not already available");
    return;
  case CANNOT_PUT_VECTOR  :
    Print(where, "Could not store vector");
    return;
  case PREPARE_NOT_OK     :
    Print(where, "DefineParameters must be called prior to this function");
    return;
  case BASIS_NOT_OK       :
    Print(where, "Arnoldi basis is not available");
    return;
  case VALUES_NOT_OK      :
    Print(where, "Eigenvalues are not available");
    return;
  case VECTORS_NOT_OK     :
    Print(where, "Eigenvectors are not available");
    return;
  case SCHUR_NOT_OK       :
    Print(where, "Schur vectors are not available");
    return;
  case RESID_NOT_OK       :
    Print(where, "Residual vector is not available");
    return;
  case MATRIX_IS_SINGULAR :
    Print(where, "Matrix is singular and could not be factored");
    return;
  case DATA_UNDEFINED     :
    Print(where, "Matrix data was not defined");
    return;
  case INSUFICIENT_MEMORY :
    Print(where, "fill-in factor must be increased");
    return;
  case NOT_SQUARE_MATRIX  :
    Print(where, "Matrix must be square to be factored");
    return;
  case NOT_FACTORED_MATRIX:
    Print(where, "Matrix must be factored before solving a system");
    return;
  case INCOMPATIBLE_SIZES :
    Print(where, "Matrix dimensions must agree");
    return;
  case DIFFERENT_TRIANGLES:
    Print(where, "A.uplo and B.uplo must be equal");
    return;
  case INCONSISTENT_DATA  :
    Print(where, "Matrix data contain inconsistencies");
    return;
  case CANNOT_READ_FILE   :
    Print(where, "Data file could not be read");
    return;
  case CANNOT_OPEN_FILE   :
    Print(where, "Invalid path or filename");
    return;
  case WRONG_MATRIX_TYPE  :
    Print(where, "Wrong matrix type");
    return;
  case WRONG_DATA_TYPE    :
    Print(where, "Wrong data type");
    return;
  case RHS_IGNORED        :
    Print(where, "RHS vector will be ignored");
    return;
  case UNEXPECTED_EOF     :
    Print(where, "Unexpected end of file");
    return;
  case NCV_OUT_OF_BOUNDS  :
    Print(where, "'ncv' is out of bounds");
    return;
  case MAXIT_NON_POSITIVE :
    Print(where, "'maxit' must be greater than zero");
    return;
  case MAX_ITERATIONS     :
    Print(where, "Maximum number of iterations taken");
    return;
  case NO_SHIFTS_APPLIED  :
    Print(where, "No shifts could be applied during a cycle of IRAM iteration");
    return;
  case CHANGING_AUTOSHIFT :
    Print(where, "Turning to automatic selection of implicit shifts");
    return;
  case DISCARDING_FACTORS :
    Print(where, "Factors L and U were not copied. Matrix must be factored");
    return;
  case GENERIC_WARNING    :
  default: ;
    Print(where, "There is something wrong");
    return;
  }

} // Set.

ArpackError::ErrorCode ArpackError::code = NO_ERRORS;
// "code" initialization.

#endif // ARERROR_H
