# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([test/decoupled/1p/test_diffusion.cc])
AM_CONFIG_HEADER([config.h])

AC_CHECK_PROGS([TEX4HT], [tex4ht], [true])
AC_CHECK_PROGS([MK4HT], [mk4ht], [true])
AC_CHECK_PROGS([T4HT], [t4ht], [true])
AM_CONDITIONAL([TEX4HT], [test "x$TEX4HT" != xtrue])
AC_CHECK_PROGS([CONVERT], [convert], [false])
AM_CONDITIONAL([CONVERT], [test "x$CONVERT" != xfalse])
AM_CONDITIONAL([BUILD_HANDBOOK], [test -a "stamp-vc"])

AC_CHECK_HEADER([valgrind/memcheck.h], 
                [HAVE_VALGRIND_H="1"],
                AC_MSG_WARN([valgrind/memcheck.h not found]))
if test x"$HAVE_VALGRIND_H" = "x1"; then
   AC_DEFINE(HAVE_VALGRIND,1,[Define whether the valgrind header files for client requests are present.])
fi

# we need no more than the standard DE-stuff
# this module depends on dune-pdelab
# this implies checking for [dune-common], [dune-grid], [dune-localfunctions], [dune-istl], [dune-pdelab]
DUNE_CHECK_ALL

# the Boost c++ template libraries
AX_BOOST_BASE([1.33.1])

AC_DEFUN([SET_PARDISO], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_pardiso_ok=no

AC_ARG_WITH(pardiso,
	[AC_HELP_STRING([--with-pardiso=<lib>], [use PARDISO library <lib>])])
case $with_pardiso in
	yes | "") ;;
	no) acx_pardiso_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) PARDISO_LIBS="$with_pardiso" ;;
	*) PARDISO_LIBS="-l$with_pardiso" ;;
esac

# Get fortran linker names of PARDISO functions to check for.
AC_F77_FUNC(pardisoinit)

acx_pardiso_save_LIBS="$LIBS"
LIBS="$BLAS_LIBS $LIBS $FLIBS"

# First, check PARDISO_LIBS environment variable
if test $acx_pardiso_ok = no; then
if test "x$PARDISO_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$PARDISO_LIBS $LIBS"
	AC_MSG_CHECKING([for $pardisoinit in $PARDISO_LIBS])
	AC_TRY_LINK_FUNC($pardisoinit, [acx_pardiso_ok=yes], [PARDISO_LIBS=""])
	AC_MSG_RESULT($acx_pardiso_ok)
	LIBS="$save_LIBS"
fi
fi

AC_SUBST(PARDISO_LIBS)

LIBS="$acx_pardiso_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_pardiso_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_PARDISO,1,[Define if you have a PARDISO library.]),[$1])
        :
else
        acx_pardiso_ok=no
        $2
fi
])dnl SET_PARDISO

SET_PARDISO

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, "$ALL_PKG_CPPFLAGS" )
AC_SUBST(AM_LDFLAGS, "$ALL_PKG_LDFLAGS" )
LIBS="$DUNE_LIBS $ALL_PKG_LIBS $PARDISO_LIBS"

AC_CONFIG_FILES([dumux.pc
    Makefile
    doc/Makefile 
    doc/doxygen/Makefile 
    doc/handbook/Makefile
    dumux/Makefile 
    dumux/boxmodels/Makefile 
    dumux/boxmodels/1p/Makefile 
    dumux/boxmodels/1p2c/Makefile 
    dumux/boxmodels/2p/Makefile 
    dumux/boxmodels/2p2c/Makefile 
    dumux/boxmodels/2p2cni/Makefile 
    dumux/boxmodels/2pni/Makefile 
    dumux/boxmodels/boxscheme/Makefile 
    dumux/boxmodels/pdelab/Makefile 
    dumux/boxmodels/richards/Makefile 
    dumux/common/Makefile 
    dumux/decoupled/Makefile 
    dumux/decoupled/2p/Makefile 
    dumux/decoupled/2p/diffusion/Makefile 
    dumux/decoupled/2p/diffusion/fv/Makefile 
    dumux/decoupled/2p/diffusion/fvmpfa/Makefile 
    dumux/decoupled/2p/diffusion/mimetic/Makefile 
    dumux/decoupled/2p/impes/Makefile 
    dumux/decoupled/2p/transport/Makefile 
    dumux/decoupled/2p/transport/fv/Makefile 
    dumux/decoupled/common/Makefile 
    dumux/io/Makefile
    dumux/material/Makefile 
    dumux/material/binarycoefficients/Makefile 
    dumux/material/components/Makefile 
    dumux/material/fluidmatrixinteractions/Makefile 
    dumux/material/fluidsystems/Makefile 
    dumux/material/spatialparameters/Makefile 
    dumux/nonlinear/Makefile 
    m4/Makefile
    test/Makefile
    test/boxmodels/Makefile 
    test/boxmodels/1p/Makefile 
    test/boxmodels/2p/Makefile 
    test/boxmodels/2pni/Makefile
    test/boxmodels/1p2c/Makefile 
    test/boxmodels/2p2c/Makefile
    test/boxmodels/2p2cni/Makefile
    test/boxmodels/richards/Makefile  
    test/common/Makefile
    test/common/pardiso/Makefile
    test/common/propertysystem/Makefile
    test/common/spline/Makefile
    test/decoupled/Makefile
    test/decoupled/1p/Makefile
    test/decoupled/2p/Makefile
    tutorial/Makefile
])

AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
echo "Pardiso..........: $acx_pardiso_ok"
echo 
