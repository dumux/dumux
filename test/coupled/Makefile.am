# tests where program to build and program to run are equal
NORMALTESTS = test_coupled test_stokesdarcy test_boxstokesdarcy

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_coupleddir = $(includedir)/dune/dumux/test/coupled
test_coupled_HEADERS = boxdiffusion.hh boxdiffusionjacobian.hh diffusionparameters.hh 
test_coupled_SOURCES = test_coupled.cc

test_stokesdarcy_SOURCES = test_stokesdarcy.cc

test_boxstokesdarcy_SOURCES = test_boxstokesdarcy.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
