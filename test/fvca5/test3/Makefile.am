# tests where program to build and program to run are equal
NORMALTESTS = fvca5_test3

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

fvca5_test3dir = $(includedir)/dune/dumux/test/fvca5/test3
fvca5_test3_HEADERS = fvca5test3problem.hh
fvca5_test3_SOURCES = fvca5_test3.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
