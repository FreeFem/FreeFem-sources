#!/usr/bin/perl
# Runs all regression tests on all compiled FreeFem++ versions
# ------------------------------------------------------------

# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 19/5/04
# $Id$

use strict;

# All existing FreeFem++ versions for this platform
my @progs;
push @progs,"$ENV{PROGLOC}/std/FreeFem++$ENV{EXEEXT}";
push @progs,"$ENV{PROGLOC}/nw/FreeFem++-nw$ENV{EXEEXT}";
push @progs,"$ENV{PROGLOC}/x11/FreeFem++-x11$ENV{EXEEXT}"
  if $ENV{X11PROG} ne "";
push @progs,"$ENV{PROGLOC}/glx/FreeFem++-glx$ENV{EXEEXT}"
  if $ENV{GLXPROG} ne "";
push @progs,"$ENV{PROGLOC}/agl/FreeFem++-agl$ENV{EXEEXT}"
  if $ENV{AGLPROG} ne "";

# Loop on all available FreeFem programs
foreach my $ff (@progs) {

  print "checking $ff\n";

  # Run FreeFem++
  system "$ff all.edp > regtests.log";
  die if $?;
}
