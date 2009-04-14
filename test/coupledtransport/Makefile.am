# tests where program to build and program to run are equal
NORMALTESTS = test_boxcoupledtransport

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_boxcoupledtransport_SOURCES = test_boxcoupledtransport.cc
AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
