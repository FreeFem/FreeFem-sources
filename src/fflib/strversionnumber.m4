#include "config-wrapper.h"
#include "strversionnumber.hpp"
#include <cstdlib>
using namespace std;
#define TOSTRING1(i) #i
#define TOSTRING(i) TOSTRING1(i)

#include <sstream>
using namespace std;

double VersionNumber(){
  return VersionFreeFempp;
}

string StrVersionNumber(){
  std::ostringstream buffer;
  buffer.precision(8);
  buffer<<VersionNumber();
  return buffer.str()+" (date VersionFreeFemDate)" ;
}
