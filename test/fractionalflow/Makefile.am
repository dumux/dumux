# tests where program to build and program to run are equal
NORMALTESTS = test_fractionalflow test_decoupled2p2c test_multiphysics2p2c

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_fractionalflow.cc test_decoupled2p2c.cc

test_fractionalflowdir = $(includedir)/dune/dumux/test/fractionalflow
test_fractionalflow_SOURCES = test_fractionalflow.cc

test_decoupled2p2c_SOURCES = test_decoupled2p2c.cc

test_multiphysics2p2c_SOURCES = test_multiphysics2p2c.cc

#AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
