# tests where program to build and program to run are equal
NORMALTESTS = test_dimltdimworld

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_dimltdimworlddir = $(includedir)/dune/dumux/test/dimltdimworld
test_dimltdimworld_SOURCES = test_dimltdimworld.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
