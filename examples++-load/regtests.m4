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
ONETEST( LapDG3)
ONETEST( LapDG4)
ONETEST( LaplaceP3)
ONETEST( LaplaceP4)
ONETEST( bilapMorley)
ONETEST( plot-fb-P3)
ONETEST( plot-fb-P3dc)
ONETEST( plot-fb-P4)
ONETEST( plot-fb-P4dc)
ONETEST( splitmesh3)
ONETEST( splitmesh6)
ONETEST( testFE-PkEdge)
ONETEST( testFE)
ONETEST( testFEMorley)


