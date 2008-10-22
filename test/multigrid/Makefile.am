dist_noinst_DATA = test_multigrid.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_multigrid

test_multigriddir = $(includedir)/dune/dumux/test/multigrid
test_multigrid_HEADERS = finiteelementerror.hh mgpreconditioner.hh testproblem.hh 
test_multigrid_SOURCES = test_multigrid.cc
test_multigrid_DEPENDENCIES = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
