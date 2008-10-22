dist_noinst_DATA = fvca5_test1.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = fvca5_test1

fvca5_test1dir = $(includedir)/dune/dumux/test/fvca5/test1
fvca5_test1_HEADERS = fvca5test1problem.hh
fvca5_test1_SOURCES = fvca5_test1.cc
fvca5_test1_DEPENDENCIES = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
