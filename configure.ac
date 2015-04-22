# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.62])
DUNE_AC_INIT # gets module version from dune.module file
AC_CONFIG_SRCDIR([dumux/common/spline.hh])
AC_CONFIG_HEADERS([config.h])

DUNE_CHECK_ALL
DUMUX_CHECKS_PRIVATE

AC_CONFIG_FILES([dumux.pc
    Makefile
    doc/Makefile
    doc/doxygen/Makefile
    doc/doxygen/Doxyfile
    doc/handbook/Makefile
    dumux/Makefile
    dumux/common/Makefile
    dumux/decoupled/Makefile
    dumux/decoupled/1p/Makefile
    dumux/decoupled/1p/diffusion/Makefile
    dumux/decoupled/1p/diffusion/fv/Makefile
    dumux/decoupled/2p/Makefile
    dumux/decoupled/2p/diffusion/Makefile
    dumux/decoupled/2p/diffusion/fv/Makefile
    dumux/decoupled/2p/diffusion/fvmpfa/Makefile
    dumux/decoupled/2p/diffusion/fvmpfa/lmethod/Makefile
    dumux/decoupled/2p/diffusion/fvmpfa/omethod/Makefile
    dumux/decoupled/2p/diffusion/mimetic/Makefile
    dumux/decoupled/2p/impes/Makefile
    dumux/decoupled/2p/transport/Makefile
    dumux/decoupled/2p/transport/fv/Makefile
    dumux/decoupled/2p2c/Makefile
    dumux/decoupled/common/Makefile
    dumux/decoupled/common/fv/Makefile
    dumux/decoupled/common/fv/mpfa/Makefile
    dumux/decoupled/common/mimetic/Makefile
    dumux/freeflow/Makefile
    dumux/freeflow/stokes/Makefile
    dumux/freeflow/stokesnc/Makefile
    dumux/freeflow/stokesncni/Makefile
    dumux/geomechanics/Makefile
    dumux/geomechanics/elastic/Makefile
    dumux/geomechanics/el1p2c/Makefile
    dumux/geomechanics/el2p/Makefile
    dumux/implicit/Makefile
    dumux/implicit/1p/Makefile
    dumux/implicit/1p2c/Makefile
    dumux/implicit/2p/Makefile
    dumux/implicit/2p2c/Makefile
    dumux/implicit/2pdfm/Makefile
    dumux/implicit/3p/Makefile
    dumux/implicit/3p3c/Makefile
    dumux/implicit/box/Makefile
    dumux/implicit/cellcentered/Makefile
    dumux/implicit/common/Makefile
    dumux/implicit/adaptive/Makefile
    dumux/implicit/co2/Makefile
    dumux/implicit/nonisothermal/Makefile
    dumux/implicit/mpnc/Makefile
    dumux/implicit/mpnc/diffusion/Makefile
    dumux/implicit/mpnc/energy/Makefile
    dumux/implicit/mpnc/mass/Makefile
    dumux/implicit/richards/Makefile
    dumux/io/Makefile
    dumux/linear/Makefile
    dumux/material/Makefile
    dumux/material/binarycoefficients/Makefile
    dumux/material/components/Makefile
    dumux/material/components/iapws/Makefile
    dumux/material/fluidmatrixinteractions/Makefile
    dumux/material/fluidmatrixinteractions/1p/Makefile
    dumux/material/fluidmatrixinteractions/2p/Makefile
    dumux/material/fluidmatrixinteractions/2pia/Makefile
    dumux/material/fluidmatrixinteractions/3p/Makefile
    dumux/material/fluidmatrixinteractions/mp/Makefile
    dumux/material/spatialparams/Makefile
    dumux/material/fluidstates/Makefile
    dumux/material/fluidsystems/Makefile
    dumux/material/constraintsolvers/Makefile
    dumux/material/eos/Makefile
    dumux/multidomain/Makefile
    dumux/multidomain/2cstokes2p2c/Makefile
    dumux/multidomain/2cnistokes2p2cni/Makefile
    dumux/multidomain/common/Makefile
    dumux/multidomain/couplinglocalresiduals/Makefile
    dumux/nonlinear/Makefile
    dumux/parallel/Makefile
    m4/Makefile
    patches/Makefile
    test/Makefile
    test/implicit/Makefile
    test/implicit/1p/Makefile
    test/implicit/1p2c/Makefile
    test/implicit/2p/Makefile
    test/implicit/2p2c/Makefile
    test/implicit/2pdfm/Makefile
    test/implicit/3p/Makefile
    test/implicit/3p3c/Makefile
    test/implicit/co2/Makefile
    test/implicit/mpnc/Makefile
    test/implicit/richards/Makefile
    test/io/Makefile
    test/io/gnuplotinterface/Makefile
    test/common/Makefile
    test/common/generalproblem/Makefile
    test/common/propertysystem/Makefile
    test/common/spline/Makefile
    test/decoupled/Makefile
    test/decoupled/1p/Makefile
    test/decoupled/2p/Makefile
    test/decoupled/2p2c/Makefile
    test/freeflow/Makefile
    test/freeflow/navierstokes/Makefile
    test/freeflow/stokes/Makefile
    test/freeflow/stokes2c/Makefile
    test/freeflow/stokes2cni/Makefile
    test/geomechanics/Makefile
    test/geomechanics/elastic/Makefile
    test/geomechanics/el1p2c/Makefile
    test/geomechanics/el2p/Makefile
    test/material/Makefile
    test/material/fluidmatrixinteractions/Makefile
    test/material/fluidmatrixinteractions/2p/Makefile
    test/material/fluidsystems/Makefile
    test/material/immiscibleflash/Makefile
    test/material/ncpflash/Makefile
    test/material/pengrobinson/Makefile
    test/material/tabulation/Makefile
    test/multidomain/Makefile
    test/multidomain/2cstokes2p2c/Makefile
    test/multidomain/2cnistokes2p2cni/Makefile
    test/references/Makefile
    tutorial/Makefile
])

# implicitly set the Dune-flags everywhere
# we need single quotes, and $(..) here because top_srcdir is only defined in the Makefile
# and not in configure
AC_SUBST(AM_CPPFLAGS, '-I$(top_srcdir) $(ALL_PKG_CPPFLAGS)' )
AC_SUBST(AM_LDFLAGS, "$ALL_PKG_LDFLAGS" )

LIBS="$DUNE_LIBS $ALL_PKG_LIBS"

AC_OUTPUT

# finally print the summary information
DUNE_SUMMARY_ALL
