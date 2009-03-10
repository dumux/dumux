# tests where program to build and program to run are equal
NORMALTESTS = test_star2dgf

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_star2dgf.cc

test_star2dgfdir = $(includedir)/dune/dumux/test/star2dgf
test_star2dgf_HEADERS = checkgeometryinfather.cc checkindexset.cc checkintersectionit.cc \
	gridcheck.cc
test_star2dgf_SOURCES = test_star2dgf.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
