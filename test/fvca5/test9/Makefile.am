# tests where program to build and program to run are equal
NORMALTESTS = fvca5_test9

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

fvca5_test9dir = $(includedir)/dune/dumux/test/fvca5/test9
fvca5_test9_HEADERS = fvca5test9problem.hh
fvca5_test9_SOURCES = fvca5_test9.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
