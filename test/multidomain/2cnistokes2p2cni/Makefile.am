if HAVE_DUNE_MULTIDOMAIN
check_PROGRAMS = test_2cnistokes2p2cni
else
headercheck_IGNORE = \
  2cnistokes2p2cniproblem.hh \
  2cnistokes2p2cnispatialparams.hh \
  2p2cnisubproblem.hh \
  stokes2cnisubproblem.hh
endif

test_2cnistokes2p2cni_SOURCES = test_2cnistokes2p2cni.cc

EXTRA_DIST = CMakeLists.txt

include $(top_srcdir)/am/global-rules
