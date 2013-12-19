#ifndef ASSERTION_HPP_
#define ASSERTION_HPP_
//  to compile all assertion
//#define ASSERTION
// to remove all the assert  
//#define NDEBUG
#ifndef ASSERTION
#define ASSERTION(i)  ((void )  0)
#else
#include <cassert>
#undef ASSERTION
#define ASSERTION(i) assert(i)
#endif
#endif
