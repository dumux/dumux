# tests where program to build and program to run are equal
NORMALTESTS = fvca5_test6

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

fvca5_test6dir = $(includedir)/dune/dumux/test/fvca5/test6
fvca5_test6_HEADERS = fvca5test6problem.hh
fvca5_test6_SOURCES = fvca5_test6.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
