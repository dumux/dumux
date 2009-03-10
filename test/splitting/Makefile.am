# tests where program to build and program to run are equal
NORMALTESTS = test_splitting

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_splitting.cc

test_splitting_SOURCES = test_splitting.cc

include $(top_srcdir)/am/global-rules
