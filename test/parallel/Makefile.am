dist_noinst_DATA = test_parallel.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_parallel

test_parallel_SOURCES = test_parallel.cc
test_parallel_CXXFLAGS = $(MPI_CPPFLAGS) 
test_parallel_LDADD = $(MPI_LDFLAGS) 

include $(top_srcdir)/am/global-rules
