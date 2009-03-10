# tests where program to build and program to run are equal
NORMALTESTS = test_pardiso

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_pardiso.cc

test_pardiso_SOURCES = test_pardiso.cc

include $(top_srcdir)/am/global-rules
