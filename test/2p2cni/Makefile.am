dist_noinst_DATA = test_2p2cni.cc test_co2.cc

noinst_PROGRAMS = test_2p2cni test_co2 test_co2_parallel 

test_2p2cnidir = $(includedir)/dune/dumux/test/2p2cni
test_2p2cni_HEADERS = brineco2problem.hh co2_soilproperties.hh waterairproblem.hh
test_2p2cni_SOURCES = test_2p2cni.cc 
test_co2_SOURCES = test_co2.cc
test_co2_parallel_SOURCES = test_co2_parallel.cc
test_co2_parallel_CXXFLAGS = $(MPI_CPPFLAGS) 
test_co2_parallel_LDADD = $(MPI_LDFLAGS) 

include $(top_srcdir)/am/global-rules
