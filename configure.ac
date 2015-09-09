#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.64)
AC_INIT([Archimedes],[2.0.1],[jeanmichel.sellier@gmail.com])
AC_CONFIG_SRCDIR([src/archimedes.c])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE([gnu dist-bzip2])

# Checks for programs.
AC_PROG_CXX

# save CFLAGS because AC_PROG_CC want to put '-g -O2' if void
save_cflags="$CFLAGS"
AC_PROG_CC
CFLAGS="$save_cflags"

AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

AC_CHECK_PROGS(LATEX, latex latex2e, no)
#AC_CHECK_PROGS(PDFLATEX, pdflatex, no)
AC_CHECK_PROGS(DVIPS, dvips, no)
AC_CHECK_PROGS(DVIPDF, dvipdf, no)


# Checks for libraries.
AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([memory.h stdlib.h string.h strings.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_STRUCT_TM

# Checks for library functions.
AC_CHECK_FUNCS([memset pow sqrt])

dnl Wether we want to create the LaTeX manual
MANUAL=nomanual
AC_ARG_ENABLE(manual,dnl
[  --enable-manual    build LaTeX manual], [ with_manual=$enableval ])

dnl only if latex/2e is found
if test "$with_manual" = yes -a "$DVIPDF" != no ; then
  MANUAL=manual
fi



AC_SUBST(LATEX)
#AC_SUBST(PDFLATEX)
AC_SUBST(MANUAL)
AC_SUBST(DVIPS)
AC_SUBST(DVIPDF)


AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 src/Makefile])
AC_OUTPUT

if test "$MANUAL" = manual; then
  AC_MSG_RESULT([  Build LaTeX manual])
fi
