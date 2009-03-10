# tests where program to build and program to run are equal
NORMALTESTS = test_multigrid

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_multigrid.cc

test_multigriddir = $(includedir)/dune/dumux/test/multigrid
test_multigrid_HEADERS = finiteelementerror.hh mgpreconditioner.hh testproblem.hh 
test_multigrid_SOURCES = test_multigrid.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
