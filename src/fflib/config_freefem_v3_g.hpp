#ifndef INIT_FREEFEM_V3_
#define INIT_FREEFEM_V3_


#ifndef VERBOSE
// #define VERBOSE
#endif


//  bamg compilation FLAG
// ----------------------
//   the freefem+ verson 
//#define NDEBUG
//#define TEST 100
//  pour bamg en version DEBUG 
//#define DEBUG
#define EIGENVALUE

#ifndef DRAWING
#define DRAWING
#endif

#ifndef LONG_LONG
#ifndef __INTEL__
#define LONG_LONG
#endif
#endif

// RNM  compilation FLAG
// ---------------------
#ifndef CHECK_KN
#define CHECK_KN
#endif

 // to use the umfpack linear solver library
#define UMFPACK


//  virtual machine exec type checking  flag (very slow)
//  pour faire un chech dynamique de tous les type du laguage
//  --------------------------------
#ifndef WITHCHECK
//#define WITHCHECK
#endif
#endif
