# tests where program to build and program to run are equal
NORMALTESTS = test_coupledcapillary test_stokesdarcy test_boxstokesdarcy test_boxstokesdarcytransport

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_coupledcapillary_SOURCES = test_coupledcapillary.cc
test_stokesdarcy_SOURCES = test_stokesdarcy.cc
test_boxstokesdarcy_SOURCES = test_boxstokesdarcy.cc
test_boxstokesdarcytransport_SOURCES = test_boxstokesdarcytransport.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
