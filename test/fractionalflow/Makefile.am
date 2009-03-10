# tests where program to build and program to run are equal
NORMALTESTS = test_fractionalflow test_decoupled2p2c test_heterogeneoussoil

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_fractionalflowdir = $(includedir)/dune/dumux/test/fractionalflow
test_fractionalflow_HEADERS = heterogeneoussoil_testproblem.hh soilproperties.hh
test_fractionalflow_SOURCES = test_fractionalflow.cc

test_decoupled2p2c_SOURCES = test_decoupled2p2c.cc

test_heterogeneoussoil_SOURCES = test_heterogeneoussoil.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
