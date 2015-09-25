# -*- Autoconf -*-
DUNE_AC_INIT
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([dumux/common/spline.hh])
AM_CONFIG_HEADER([config.h])

DUMUX_CHECK_ALL

AC_CONFIG_FILES([dumux.pc
    Makefile
    appl/Makefile
    appl/lecture/Makefile
    appl/lecture/msm/Makefile
    appl/lecture/msm/1p2cvs2p/Makefile
    appl/lecture/msm/buckleyleverett/Makefile
    appl/lecture/msm/mcwhorter/Makefile
    doc/Makefile 
    doc/doxygen/Makefile 
    doc/doxygen/Doxyfile
    doc/handbook/Makefile
    dumux/Makefile 
    dumux/boxmodels/Makefile 
    dumux/boxmodels/1p/Makefile 
    dumux/boxmodels/1p2c/Makefile 
    dumux/boxmodels/2p/Makefile 
    dumux/boxmodels/2p2c/Makefile 
    dumux/boxmodels/2p2cni/Makefile 
    dumux/boxmodels/2pni/Makefile 
    dumux/boxmodels/common/Makefile 
    dumux/boxmodels/richards/Makefile 
    dumux/common/Makefile 
    dumux/decoupled/Makefile
    dumux/decoupled/1p/Makefile
    dumux/decoupled/1p/diffusion/Makefile
    dumux/decoupled/1p/diffusion/fv/Makefile 
    dumux/decoupled/2p/Makefile 
    dumux/decoupled/2p/diffusion/Makefile 
    dumux/decoupled/2p/diffusion/fv/Makefile 
    dumux/decoupled/2p/diffusion/fvmpfa/Makefile 
    dumux/decoupled/2p/diffusion/mimetic/Makefile 
    dumux/decoupled/2p/impes/Makefile 
    dumux/decoupled/2p/transport/Makefile 
    dumux/decoupled/2p/transport/fv/Makefile
    dumux/decoupled/2p2c/Makefile 
    dumux/decoupled/common/Makefile 
    dumux/io/Makefile
    dumux/material/Makefile 
    dumux/material/binarycoefficients/Makefile 
    dumux/material/components/Makefile 
    dumux/material/components/iapws/Makefile 
    dumux/material/fluidmatrixinteractions/Makefile 
    dumux/material/fluidmatrixinteractions/2p/Makefile 
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
    test/decoupled/2p2c/Makefile
    test/material/Makefile
    test/material/tabulation/Makefile
    tutorial/Makefile
])

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, "$ALL_PKG_CPPFLAGS" )
AC_SUBST(AM_LDFLAGS, "$ALL_PKG_LDFLAGS" )

LIBS="$DUNE_LIBS $ALL_PKG_LIBS"

AC_OUTPUT

# finally print the summary information
DUMUX_SUMMARY_ALL
