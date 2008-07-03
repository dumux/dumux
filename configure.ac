# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([test/diffusion/test_diffusion.cc])
AM_CONFIG_HEADER([config.h])
#AC_SUBST([LOCAL_LIBS], '$(top_srcdir)/lib/libdunemux.la')
AC_SUBST([LOCAL_LIBS], '$(top_srcdir)/dumux/libdumux.la')

# we need no more than the standard DE-stuff
# this module depends on dune-common dune-disc dune-fem dune-grid dune-istl
# this implies checking for [dune-common], [dune-istl], [dune-grid], [dune-disc], [dune-fem]
DUNE_CHECK_ALL([dune-common], [dune-istl], [dune-grid], [dune-disc])

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
AC_SUBST(AM_CPPFLAGS, $ALL_PKG_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $ALL_PKG_LDFLAGS)
LIBS="$DUNE_LIBS $ALL_PKG_LIBS $PARDISO_LIBS"

AC_CONFIG_FILES([Makefile
    doc/Makefile 
    dumux/Makefile 
    test/Makefile
    doc/doxygen/Makefile 
    dumux/diffusion/Makefile 
    dumux/functions/Makefile 
    dumux/io/Makefile 
    dumux/material/Makefile 
    dumux/onedinndgrid/Makefile
    dumux/operators/Makefile 
    dumux/shapefunctions/Makefile 
    dumux/timedisc/Makefile 
    dumux/transport/Makefile 
    dumux/fractionalflow/Makefile 
    dumux/twophase/Makefile 
    dumux/2penergy/Makefile
    dumux/parabolic/Makefile 
    dumux/nonlinear/Makefile 
    test/boundaryid/Makefile
    test/box3d/Makefile
    test/brinkman/Makefile
    test/capillary/Makefile
    test/co2/Makefile
    test/diffusion/Makefile
    test/dimltdimworld/Makefile
    test/elementid/Makefile
    test/fractionalflow/Makefile
    test/fvca5/Makefile
    test/gravity/Makefile
    test/multigrid/Makefile
    test/newton/Makefile
    test/parallel/Makefile
    test/pardiso/Makefile
    test/splitting/Makefile
    test/star2dgf/Makefile
    test/star2uggrid/Makefile
    test/transport/Makefile
    test/twophase/Makefile
    test/upscaledvelocity/Makefile
    test/vtk/Makefile
    dumux/diffusion/fv/Makefile 
    dumux/diffusion/fe/Makefile 
    dumux/diffusion/mimetic/Makefile 
    dumux/diffusion/problems/Makefile 
    dumux/parabolic/fe/Makefile 
    dumux/parabolic/fv/Makefile 
    dumux/parabolic/problems/Makefile 
    dumux/twophase/fv/Makefile 
    dumux/twophase/problems/Makefile 
    dumux/2penergy/fv/Makefile
    dumux/transport/fv/Makefile 
    dumux/transport/problems/Makefile 
    dumux/fractionalflow/impes/Makefile 
    dumux/shapefunctions/cr/Makefile
    test/fvca5/test1/Makefile
    test/fvca5/test2/Makefile
    test/fvca5/test3/Makefile
    test/fvca5/test4/Makefile
    test/fvca5/test5/Makefile
    test/fvca5/test6/Makefile
    test/fvca5/test7/Makefile
    test/fvca5/test8/Makefile
    test/fvca5/test9/Makefile
    benchmarkproblems/Makefile
    benchmarkproblems/fractionalflow/Makefile
    benchmarkproblems/fullycoupled/Makefile
    benchmarkproblems/fullycoupled/2p/Makefile
    benchmarkproblems/fractionalflow/2p/Makefile
    benchmarkproblems/fullycoupled/2p/buckleyleverett/Makefile
    benchmarkproblems/fractionalflow/2p/buckleyleverett/Makefile
    benchmarkproblems/fullycoupled/2p/fivespot/Makefile
    benchmarkproblems/fractionalflow/2p/fivespot/Makefile    
    benchmarkproblems/fullycoupled/2p/convectiondiffusion/Makefile  
    benchmarkproblems/fractionalflow/2p/convectiondiffusion/Makefile      
    benchmarkproblems/fullycoupled/2p/mcwhorter/Makefile
    benchmarkproblems/fractionalflow/2p/mcwhorter/Makefile
    dune-mux.pc])

AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
echo "Pardiso..........: $acx_pardiso_ok"
echo 
