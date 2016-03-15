# Checks whether a compiler accepts a given flag
# ----------------------------------------------

# $1 = compiler name
# $2 = flag
# $3 = make macro containing flags for that compiler
# $4 =  exec is true ..
# Note: changes AC_LANG()

AC_DEFUN([CHECK_COMPILE_FLAG],
	[AC_MSG_CHECKING(whether the $1 compiler accepts $2)
	check_save_flags="$$3"
	AC_LANG_PUSH($1)
	$3="$$3 $2"

	# The program needs to contain something for the test source
	# file to be created by autoconf.

	# Some options really need to be linked (not only compiled) to
	# check whether they work.

	AC_LINK_IFELSE([ifelse($1,Fortran 77,
[       program x
       end],
			[AC_LANG_PROGRAM])],
		check_flag_ok=yes,
		check_flag_ok=no)
	AC_MSG_RESULT($check_flag_ok)
	if test "$check_flag_ok" = no;
	then
		$3="$check_save_flags"
	fi
	if test -n "$4" ;
         then 
	        $4="$check_flag_ok"
	fi

	AC_LANG_POP($1)
])
