if HAVE_DUNE_MULTIDOMAIN
check_PROGRAMS = test_2cstokes2p2c
else
headercheck_IGNORE = \
  2cstokes2p2cproblem.hh \
  2cstokes2p2cspatialparams.hh \
  2p2csubproblem.hh \
  stokes2csubproblem.hh
endif

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST := $(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_2cstokes2p2c_SOURCES = test_2cstokes2p2c.cc

include $(top_srcdir)/am/global-rules
