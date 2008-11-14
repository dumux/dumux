dist_noinst_DATA = test_parallel.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_parallel test_transport_parallel

test_paralleldir = $(includedir)/dune/dumux/test/parallel
test_parallel_HEADERS = diffusionparameters.hh parallelboxdiffusion.hh \
	parallelboxdiffusionjacobian.hh schwarzcommunication.hh schwarzpreconditioner.hh \
	simpleproblem.hh
test_parallel_SOURCES = test_parallel.cc
test_parallel_CXXFLAGS = $(MPI_CPPFLAGS) 
test_parallel_LDADD = $(MPI_LDFLAGS) 

test_transport_parallel_SOURCES = test_transport_parallel.cc
test_transport_parallel_CXXFLAGS = $(MPI_CPPFLAGS) 
test_transport_parallel_LDADD = $(MPI_LDFLAGS) 

include $(top_srcdir)/am/global-rules
