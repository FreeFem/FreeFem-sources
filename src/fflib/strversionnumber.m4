#include <config.h>
#include "strversionnumber.hpp"
#include <sstream>
using namespace std;

#define xstrg(s) strg(s)
#define strg(s) #s

double VersionNumber() {
  return VersionFreeFem;
}

string StrVersionNumber() {
  ostringstream version;
  version.precision(8);
  version << xstrg(VersionFreeFem)
          << " (VersionFreeFemDate - git GitVersion)";
  return version.str();
}
