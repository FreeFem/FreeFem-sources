#include "versionnumber.hpp"
#include <cstdlib>
using namespace std;
#define TOSTRING1(i) #i
#define TOSTRING(i) TOSTRING1(i)
const char * StrVersionNumber()
{
  
  return   TOSTRING( VersionFreeFempp ) " (date " VersionFreeFemDate ")" ;
  
}
double  VersionNumber()
{
  
  return   VersionFreeFempp ;
  
}
