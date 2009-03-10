# tests where program to build and program to run are equal
NORMALTESTS = test_minc

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_minc.cc 

test_mincdir = $(includedir)/dune/dumux/test/minc
test_minc_HEADERS = minclensproblem.hh
test_minc_SOURCES = test_minc.cc

include $(top_srcdir)/am/global-rules
