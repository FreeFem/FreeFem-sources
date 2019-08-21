#include <config.h>
#include "strversionnumber.hpp"
#include <sstream>
using namespace std;

double VersionNumber() {
  return VersionFreeFem;
}

string StrVersionNumber() {
  ostringstream version;
  version.precision(8);
  version << VersionNumber()
          << " (VersionFreeFemDate - git GitVersion)";
  return version.str();
}
