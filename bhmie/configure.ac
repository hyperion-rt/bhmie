AC_INIT(bhmie,version-dev)
AC_PROG_FC([ifort gfortran f95 g95 xlf])

if test "$FC" == f95
then
	AC_SUBST(posix_module,posix_nag)
	AC_SUBST(extra,'-Isrc/modules -mdir src/modules -f2003 -colour -ieee=full')
elif test "$FC" == ifort
then
	AC_SUBST(posix_module,posix_empty)
	AC_SUBST(extra,'-module src/modules -heap-arrays')
elif test "$FC" == gfortran
then
	AC_SUBST(posix_module,posix_empty)
	AC_SUBST(extra,'-Jsrc/modules -Isrc/modules -ffree-line-length-none')
fi

AC_OUTPUT(Makefile)
