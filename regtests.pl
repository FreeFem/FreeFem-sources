#!/usr/bin/perl
# Runs all regression tests on all compiled FreeFem++ versions
# ------------------------------------------------------------

# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 19/5/04
# $Id$

use strict;

# All existing FreeFem++ versions for this platform
my @progs;

# At present, we do not test the standard version because its window
# is too invasive (at least under X11).

#push @progs,"$ENV{PROGLOC}/std/FreeFem++$ENV{EXEEXT}";

push @progs,"$ENV{PROGLOC}/nw/FreeFem++-nw$ENV{EXEEXT}";
push @progs,"$ENV{PROGLOC}/x11/FreeFem++-x11$ENV{EXEEXT}"
  if $ENV{X11PROG} ne "";
push @progs,"$ENV{PROGLOC}/glx/FreeFem++-glx$ENV{EXEEXT}"
  if $ENV{GLXPROG} ne "";
push @progs,"$ENV{PROGLOC}/agl/FreeFem++-agl$ENV{EXEEXT}"
  if $ENV{AGLPROG} ne "";
push @progs,"$ENV{PROGLOC}/mpi/FreeFem++-mpi$ENV{EXEEXT}"
  if $ENV{MPIPROG} ne "";

# Loop on all available FreeFem programs
my $count=0;
foreach my $ff (@progs) {
  $count++;

  # Running FreeFem++ on regtests.edp (specific to regression tests),
  # otherwise on all.edp.
  my $test="all.edp";
  $test="regtests.edp" if -f "regtests.edp";
  my $cmd="$ff $test > regtests-$count.log";
  print "$cmd\n";
  system $cmd;
  die if $?;
}
