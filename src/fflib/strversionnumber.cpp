#include "strversionnumber.hpp" // [[file:strversionnumber.hpp]]
#include <cstdlib>
using namespace std;
#define TOSTRING1(i) #i
#define TOSTRING(i) TOSTRING1(i)

//#include <sstream>
#include <cstdio>
using namespace std;

double VersionNumber(){
  return VersionFreeFempp;
}

string StrVersionNumber(){
//  std::ostringstream buffer;
//  buffer.precision(8);
//  buffer<<VersionNumber();
  static char buffer[100];
  sprintf(buffer," %9f (date Mar 27 sep 2016 19:24:14 CEST)",VersionNumber());
  return buffer; //.str()+" (date Mar 27 sep 2016 19:24:14 CEST)" ;
}
