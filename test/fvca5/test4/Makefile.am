# tests where program to build and program to run are equal
NORMALTESTS = fvca5_test4

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

fvca5_test4dir = $(includedir)/dune/dumux/test/fvca5/test4
fvca5_test4_HEADERS = fvca5test4problem.hh
fvca5_test4_SOURCES = fvca5_test4.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
