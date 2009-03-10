# tests where program to build and program to run are equal
NORMALTESTS = test_newton

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_newton.cc

test_newton_SOURCES = test_newton.cc

include $(top_srcdir)/am/global-rules
