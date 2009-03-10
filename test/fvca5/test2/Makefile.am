# tests where program to build and program to run are equal
NORMALTESTS = fvca5_test2

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

fvca5_test2dir = $(includedir)/dune/dumux/test/fvca5/test2
fvca5_test2_HEADERS = fvca5test2problem.hh
fvca5_test2_SOURCES = fvca5_test2.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
