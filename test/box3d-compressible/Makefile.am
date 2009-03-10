# tests where program to build and program to run are equal
NORMALTESTS = test_box3d_compressible

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_box3d_compressibledir = $(includedir)/dune/dumux/test/box3d-compressible
test_box3d_compressible_HEADERS = boxdiffusion.hh boxdiffusionjacobian.hh diffusionparameters.hh 
test_box3d_compressible_SOURCES = test_box3d_compressible.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
