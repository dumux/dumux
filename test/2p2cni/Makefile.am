# tests where program to build and program to run are equal
NORMALTESTS = new_leakywell new_test_2p2cni test_2p2cni

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_2p2cnidir = $(includedir)/dune/dumux/test/2p2cni
test_2p2cni_HEADERS = waterairproblem.hh

new_leakywell_SOURCES = new_leakywell.cc 

new_test_2p2cni_SOURCES = new_test_2p2cni.cc 

test_2p2cni_SOURCES = test_2p2cni.cc 

include $(top_srcdir)/am/global-rules
