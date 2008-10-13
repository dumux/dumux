dist_noinst_DATA = test_star2uggrid.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_star2uggrid
test_star2uggriddir = $(includedir)/dune/dumux/test/star2uggrid
test_star2uggrid_HEADERS = checkgeometryinfather.cc checkindexset.cc checkintersectionit.cc \
	gridcheck.cc
test_star2uggrid_SOURCES = test_star2uggrid.cc
test_star2uggrid_DEPENDENCIES = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
