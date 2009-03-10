# tests where program to build and program to run are equal
NORMALTESTS = test_brinkman

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_brinkmandir = $(includedir)/dune/dumux/test/brinkman
test_brinkman_HEADERS = brinkmantestproblem.hh
test_brinkman_SOURCES = test_brinkman.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
