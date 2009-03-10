# tests where program to build and program to run are equal
NORMALTESTS = test_vtk

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_vtk.cc

test_vtk_SOURCES = test_vtk.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
