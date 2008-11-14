dist_noinst_DATA = test_star2dgf.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_star2dgf

test_star2dgfdir = $(includedir)/dune/dumux/test/star2dgf
test_star2dgf_HEADERS = checkgeometryinfather.cc checkindexset.cc checkintersectionit.cc \
	gridcheck.cc
test_star2dgf_SOURCES = test_star2dgf.cc

include $(top_srcdir)/am/global-rules
