#ifndef THROWASSERT
#define THROWASSERT
#include <iostream>

//#ifdef __INTEL__
#define cerr cout 
//#endif

#include "error.hpp"

#ifdef NDEBUG
#define throwassert(i) ((void) 0)
#else
#define throwassert(condition)  ((condition) ? ((void) 0) : throw(ErrorAssert(#condition,__FILE__, __LINE__)))
 
#undef assert
#define assert(condition) throwassert(condition)
#endif

#define InternalError(message) throw(ErrorInternal(message,__LINE__,__FILE__))
#endif
