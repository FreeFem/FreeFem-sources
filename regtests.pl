#!/usr/bin/perl
# Runs all regression tests on all compiled FreeFem++ versions
# ------------------------------------------------------------

# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 19/5/04
# $Id$

use strict;

# All FreeFem++ versions
my @progs;
push @progs,"$ENV{PROGLOC}/std/FreeFem++$ENV{EXEEXT}"
  if -e "$ENV{PROGLOC}/std/FreeFem++$ENV{EXEEXT}";
push @progs,"$ENV{PROGLOC}/nw/FreeFem++-nw$ENV{EXEEXT}"
  if -e "$ENV{PROGLOC}/nw/FreeFem++-nw$ENV{EXEEXT}";
push @progs,"$ENV{PROGLOC}/x11/FreeFem++-x11$ENV{EXEEXT}"
  if -e "$ENV{PROGLOC}/x11/FreeFem++-x11$ENV{EXEEXT}";
push @progs,"$ENV{PROGLOC}/glx/FreeFem++-glx$ENV{EXEEXT}"
  if -e "$ENV{PROGLOC}/glx/FreeFem++-glx$ENV{EXEEXT}";
push @progs,"$ENV{PROGLOC}/agl/FreeFem++-agl$ENV{EXEEXT}"
  if -e "$ENV{PROGLOC}/agl/FreeFem++-agl$ENV{EXEEXT}";
push @progs,"$ENV{PROGLOC}/pc/FreeFem++-pc$ENV{EXEEXT}"
  if -e "$ENV{PROGLOC}/pc/FreeFem++-pc$ENV{EXEEXT}";

# Loop on all available FreeFem programs
foreach my $ff (@progs) {

  print "checking $ff\n";

  # Run FreeFem++
  system "$ff all.edp > regtests.log";
  die if $?;
}
