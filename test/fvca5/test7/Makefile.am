# tests where program to build and program to run are equal
NORMALTESTS = fvca5_test7

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

fvca5_test7dir = $(includedir)/dune/dumux/test/fvca5/test7
fvca5_test7_HEADERS = fvca5test7problem.hh
fvca5_test7_SOURCES = fvca5_test7.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
