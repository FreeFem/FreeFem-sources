#!/usr/bin/perl
# Runs all regression tests on all FreeFem++ versions
# ---------------------------------------------------

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

# All tests found locally
my @tests = <*.edp>;

# Loop on all .edp files
foreach my $edp (@tests) {

# Avoid re-using "regtests.edp", which is produced locally.
  if ($edp ne "regtests.edp") {

    # Only tries the edp if it contains a value to check against
    my $file = `cat $edp`;
    if ($file =~ /Regression Test: \/(.*)\//m) {
      my $check = $1;

      # Creating temporary file with special test-only options
      open EDP,"> regtests.edp";
      print EDP "NoUseOfWait=true;\n";
      print EDP $file;
      close EDP;

      # Loop on all FreeFems
      foreach my $ff (@progs) {

	print "checking $ff on $edp\n";

	# Run FreeFem++
	my $output = `$ff regtests.edp`;
	die "$ff $edp could not be run" if $?;

	# Checks the result
	if ($output !~ /$check/) {
	  print $output;
	  die "\"$ff $edp\" does not match /$check/";
	}
      }
    }
  }
}
