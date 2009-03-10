# tests where program to build and program to run are equal
NORMALTESTS = fvca5_test8

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

fvca5_test8dir = $(includedir)/dune/dumux/test/fvca5/test8
fvca5_test8_HEADERS = fvca5test8problem.hh
fvca5_test8_SOURCES = fvca5_test8.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
