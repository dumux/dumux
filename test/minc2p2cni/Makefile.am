# tests where program to build and program to run are equal
NORMALTESTS = test_mincco2

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_mincco2.cc

test_minc2p2cnidir = $(includedir)/dune/dumux/test/minc2p2cni
test_minc2p2cni_HEADERS = brineco2problem.hh co2_soilproperties.hh waterairproblem.hh
test_mincco2_SOURCES = test_mincco2.cc

include $(top_srcdir)/am/global-rules
