# Checks whether a compiler accepts a given flag
# ----------------------------------------------

# $1 = compiler name
# $2 = flag
# $3 = make macro containing flags for that compiler

# Note: changes AC_LANG()

AC_DEFUN(CHECK_COMPILE_FLAG,
	[AC_MSG_CHECKING(whether the $1 compiler accepts $2)
	check_save_flags="$$3"
	AC_LANG_PUSH($1)
	$3="$$3 $2"

	# The program needs to contain something for the test source
	# file to be created by autoconf.

	AC_COMPILE_IFELSE([
		],
		check_flag_ok=yes,
		check_flag_ok=no)
	AC_MSG_RESULT($check_flag_ok)
	if test "$check_flag_ok" = no;
	then
		$3="$check_save_flags"
	fi

	AC_LANG_POP($1)
])
