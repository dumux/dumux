# tests where program to build and program to run are equal
NORMALTESTS = test_star2uggrid test_lgm2uggrid

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_star2uggrid.cc

test_star2uggriddir = $(includedir)/dune/dumux/test/star2uggrid
test_star2uggrid_HEADERS = checkgeometryinfather.cc checkindexset.cc checkintersectionit.cc \
	gridcheck.cc
test_star2uggrid_SOURCES = test_star2uggrid.cc

test_lgm2uggrid_SOURCES = test_lgm2uggrid.cc
test_lgm2uggrid_DEPENDENCIES = $(LOCAL_LIBS)

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
