#! /usr/bin/perl -w
# fort77 (compiler driver) script for f2c
# For use with gcc under Linux
# This code is in the public domain; use at your own risk.
# Parse options

$version = "1.14a";
$nnflag = '-Nn802';
$tmpdir = $ENV{'TMPDIR'} || '/tmp';
$cpp = 0;
$fast_math = 1;
$debug = 0;
$debugcmd = "";
push(@includes, "@INC@");
# Loop over all options; pull all options from @ARGV and put all
# arguments into @argv.	 This is needed because, apparently, UNIX
# compilers acceppt options anywhere on the command line.

while ($_ = $ARGV[0]) {
    shift;

    if (!/^-/) {
	if (/\.P$/) {
	    push(@pfiles, $_);
	}
	else {
	    push(@argv, $_);
	}
	next;
    }
    # First, the f2c options.

    if (/^-[CUuaEhRrz]$/ || /^-I[24]$/ || /^-onetrip$/ || /^-![clPR]$/ ||
	/^-ext$/ || /^-!bs$/ || /^-W[1-9][0-9]*$/ || /^-w8$/ || /^-w66$/ ||
	/^-r8$/ || /^-N[^n][0-9]+$/) {
	push (@fopts, $_);
    }
    elsif (/^-Nn[0-9]+$/) {
	$nnflag = $_;
    }

# Prototype flags for f2c

    elsif (/^-Ps?/) {
	$extract_prototypes ++;
	push (@fopts, $_);
    }

# Does somebody want to run the preprocessor?

    elsif (/^-cpp$/) {
	$cpp++;
    }

# These are common to both f2c and gcc
    elsif (/^-w$/) {
	push(@fopts, $_);
	push(@copts, $_);
    }

# This is for the linker, too...
    elsif (/^-g$/) {
	push(@fopts, $_);
	push(@copts, $_);
	push(@lopts, $_);
	$debug ++;
    }

# Special options for the different subprocesses: f for f2c step,
# p for (separate) preprocessing, c for C compiler, l for linker.
# a is also passed to the C compiler.

    elsif (/^-Wf,/) {
	push(@fopts, &parsewx($_));
    }
    elsif (/-Wp,/) {
	push(@cppopts, &parsewx($_));
    }
    elsif (/-W[ca],/) {
	push(@copts, &parsewx($_));
    }
    elsif (/-Wl,/) {
	push(@lopts,&parsewx($_));
    }

# gcc only options

# too many -f and -W options to list them all...

# First, let's see wether somebody wants to adhere to the C standard
# in Fortran.

    elsif (/^-fnofast-math$/) {
	$fast_math = 0;
    }

# The '-f' option to f2c...

    elsif (/^-f$/) {
	push(@fopts, $_);
    }
    elsif (/^-[fWUAm]/ || /^-[Ex]$/ || /^-pipe$/ ) {
	push(@copts, $_);
    }

# Includes and outputs...

    elsif (/^-I$/) {
	(@ARGV > 0) || die "$0: Missing argument to \"$_\"\n";
	push(@includes, "-I".shift);
    }
    elsif (/^-I./) {
	push(@includes, $_);
    }
    elsif (/^-o$/) {
	(@ARGV > 0) || die "$0: Missing argument to \"$_\"\n";
	$output = shift;
    }
    elsif (/^-o(.*)/) {
	$output = $1;
    }

# Optimization
    elsif (/^-O/) {
	push(@copts, $_);
	push(@lopts, $_);
	$optimize ++;
    }

# Options for both C compiler and linker

    elsif (/^-[Og]/ || /^-p$/ || /^-pg$/) {
	push(@copts, $_);
	push(@lopts, $_);
    }
    elsif (/^-[bV]$/ ) {
	(@ARGV > 0) || die "$0 : Missing argument to \"$_\"\n";
	$arg = shift;
	push(@copts, $_, $arg);
	push(@lopts, $_, $arg);
    }
    elsif (/^-[bV]./ ) {
	push(@copts, $_);
	push(@lopts, $_);
    }

# Linker only options

    elsif (/^-[lL]$/) {
	push(@lopts, $_);
	(@ARGV > 0) || die "$0: Missing argument to \"$_\"\n";
	$_ = shift;
	push(@lopts, $_);
    }
    elsif (/^-[lL]./ || /^-nostartfiles$/ || /^-static$/ || /^-shared$/ ||
	   /^-symbolic$/) {
	push(@lopts, $_);
    }
    elsif (/^-[cS]$/) {
	$compile_only = $_;
    }
    elsif (/^-D/) {
	push(@cppopts, $_);
    }
#   Are we verbose?

    elsif (/^-v$/) {
	$verbose ++;
    }

# Does somebody want to keep the C files around?

    elsif (/^-k$/) {
	$keep_c ++;
    }

    else {
	die "$0: Illegal option: $_\n";
    }

}

push(@fopts,$nnflag);
push(@copts,'-ffast-math') if $optimize && $fast_math;
push(@cppopts,@includes);
push(@fopts,@includes,"-I.");
push(@fopts, @pfiles);

if ($verbose) {
    print STDERR "$0: fort77 Version $version\n";
    if ($verbose > 1) {
	push(@copts,"-v");
	push(@lopts,"-v");
	push(@cppopts,"-v");
    }
}


@ARGV = @argv;

if ($compile_only && $output && (@ARGV>1)) {
    warn "$0: Warning: $compile_only and -o with mutiple files, ignoring -o\n";
    $output = "";
}

die "$0: No input files specified\n" unless @ARGV;

while ($_ = $ARGV[0]) {
    shift;
    $ffile = "";
    $cfile = "";
    $lfile = "";
    $basefile = "";

    if (/\.[fF]$/) {
	$ffile = $_;
	$basefile = $ffile;
    }
    elsif (/\.[cCisSm]$/ || /\.cc$/ || /\.cxx$/) {
	$cfile = $_;
	$basefile = $_;
    }
    else {
	push(@lfiles, $_);
    }
    if ($ffile) {
	&check_file_read($ffile);
	if ($keep_c) {
	    $cfile = ($ffile =~ /([^\/]*\.).$/)[0] . "c";
	}
	else {
	    $seq ++;
	    $cfile = "$tmpdir/fort77-$$-$seq.c";
	}
	if ($debug) {
	    $debugcmd = ' | /usr/bin/perl -p -e \'s/^(#line.*)""/$1"'
		. $ffile . '"/\' '
	}
	if ($cpp || ($ffile =~ /\.F$/)) {
#         Backslashes at the end of comment lines confuse cpp...
	    $pipe =  "| /lib/cpp -traditional " . 
		join(' ',@cppopts) . " | @f2c@ " .
		    join(' ',@fopts) . $debugcmd . " > $cfile";
	    print STDERR "$0: Running \"$pipe\"" if $verbose;
	    open(F2C,$pipe);

	    open (FFILE, "$ffile") || die ("$0: Cannot open $ffile: $_\n");
	    while (<FFILE>) {
		s/([cC*].*)\\$/$1/;
		print F2C $_;
	    }
	    close(FFILE);
	    close(F2C);
	    $retcode = $? / 256;

	}
	else {
	    $retcode = &mysystem("@f2c@ ".
				 join (" ",@fopts). " < ". $ffile .
				 $debugcmd . " > $cfile")/256;
	}
	if ($retcode && !$keep_c) {
	    print STDERR "$0: unlinking $cfile\n" if $verbose;
	    unlink $cfile;
	    die "$0: aborting compilation\n";
	}

# Separate the prototypes out from the C files.

	if ($extract_prototypes) {
	    $pfile = ($basefile =~ /([^\/]*\.).$/)[0] . "P";
            open(CFILE, "$cfile") || die ("$0: Cannot open $cfile\n");
# *wdh*	    while (($line = <CFILE>) &&
	    while (defined($line = <CFILE>) &&
		   ($line !~ '#ifdef P_R_O_T_O_T_Y_P_E_S\n')) {
		print $line;
	    }
	    if ($_) {
        	open(PFILE, ">$pfile") || die ("$0: Cannot open $pfile\n");
# *wdh*		while (($line = <CFILE>) && ($line !~ '#endif')) {
		while (defined($line = <CFILE>) && ($line !~ '#endif')) {
		    print PFILE $line;
		}
		close(PFILE);
	    }
	    close(CFILE);
	}
    }

# C compilation step.

    if ($cfile) {
# *wdh*	@command = ("cc",@cppopts,@copts);
	@command = ("@CC@",@cppopts,@copts);
	if ($compile_only && $output) {
	    push(@command,'-o',$output,$compile_only);
	}
	elsif ((!$compile_only) || ($compile_only eq '-c')) {
	    $lfile = ($basefile =~ /([^\/]*\.).$/)[0] . "o";
	    push(@command, '-c', '-o', $lfile);
	}
	elsif ($compile_only eq '-S') {
	    $sfile = ($basefile =~ /([^\/]*\.).$/)[0] . "s";
	    push(@command, '-S', '-o', $sfile);
	}

	push(@command,$cfile);
	$retcode = &mysystem(@command)/256;

	if ($retcode) {
	    die "$0: aborting compilation\n";
	}
	if ($ffile && !$keep_c) {
	    print STDERR "$0: unlinking $cfile\n" if $verbose;
	    unlink $cfile;
	}
	if ($lfile) {
	    push (@gener_lfiles, $lfile); push(@lfiles, $lfile);
	    $lfile = "";
	}
    }
    push (@lfiles, $lfile) if $lfile;
}


exit if $compile_only;

push (@output, "-o", $output) if $output;

$retcode = &mysystem("@CC@", @output, @lfiles, @lopts, "@LLIBDIR@","-lf2c", "-lm" );
if (@gener_lfiles) {
    print STDERR "$0: unlinking ",join(',',@gener_lfiles),"\n" if $verbose;
    unlink (@gener_lfiles);
}
exit $retcode;

# Basically a system call, except that we want to be verbose if
# necessary.

sub mysystem
{
    local (@args) = @_;
    if (@args == 1) {
	print STDERR "$0: Running \"$args[0]\"\n" if $verbose;
	system($args[0]);
    }
    else {
	print STDERR "$0: Running \"",join(' ',@args),"\"\n" if $verbose;
	system(@args);
    }
}

sub parsewx
{
    local ($str) = @_;
    local(@tmp) = split(/,/,$str);
    shift(@tmp);
    return @tmp;
}

sub check_file_read
{
    local ($name) = @_;
    open (TESTFILE,"$name") || die "Cannot open $name: $!\n";
    close(TESTFILE);
}
