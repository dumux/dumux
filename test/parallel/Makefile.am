dist_noinst_DATA = test_parallel.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_parallel

test_parallel_SOURCES = test_parallel.cc
test_parallel_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
test_parallel_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS) 
test_parallel_DEPENDENCIES = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
