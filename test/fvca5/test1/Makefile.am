# tests where program to build and program to run are equal
NORMALTESTS = fvca5_test1

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

fvca5_test1dir = $(includedir)/dune/dumux/test/fvca5/test1
fvca5_test1_HEADERS = fvca5test1problem.hh
fvca5_test1_SOURCES = fvca5_test1.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
