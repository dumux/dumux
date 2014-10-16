if HAVE_DUNE_MULTIDOMAIN
check_PROGRAMS = test_2cnistokes2p2cni
else
headercheck_IGNORE = \
  2cnistokes2p2cniproblem.hh \
  2cnistokes2p2cnispatialparams.hh \
  2p2cnisubproblem.hh \
  stokes2cnisubproblem.hh
endif

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST := $(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_2cnistokes2p2cni_SOURCES = test_2cnistokes2p2cni.cc

include $(top_srcdir)/am/global-rules
