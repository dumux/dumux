# tests where program to build and program to run are equal
NORMALTESTS = new_test_richards richards

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = new_test_richards.cc richards.cc

richards_SOURCES = richards.cc

new_test_richards_SOURCES = new_test_richards.cc

include $(top_srcdir)/am/global-rules
