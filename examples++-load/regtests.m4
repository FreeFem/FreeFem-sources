// Regression tests
// ----------------

// $Id$

include(../regtests.m4)

// The values tested here may not have a physical or mathematical
// meaning. Their main property is to gather numerical values from the
// whole domain, to be checked for consistency with previous runs.
ONETEST(load,uh[].max,0.0001)
ONETEST(testFE)
ONETEST(testFEMorley)
ONETEST(funcTemplate)

