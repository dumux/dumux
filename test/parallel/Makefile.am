# tests where program to build and program to run are equal
NORMALTESTS = test_parallel test_transport_parallel

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_parallel.cc

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

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
