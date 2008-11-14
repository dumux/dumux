dist_noinst_DATA = test_boundaryid.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_boundaryid

test_boundaryiddir = $(includedir)/dune/dumux/test/boundaryid
test_boundaryid_HEADERS = lensproblemwithid.hh checkgeometryinfather.cc checkindexset.cc \
	checkintersectionit.cc gridcheck.cc
test_boundaryid_SOURCES = test_boundaryid.cc

include $(top_srcdir)/am/global-rules
