# tests where program to build and program to run are equal
check_PROGRAMS = test_diffusion test_diffusion3d test_1p

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_diffusion_SOURCES = test_diffusion.cc

test_diffusion3d_SOURCES = test_diffusion3d.cc

test_1p_SOURCES = test_1p.cc

include $(top_srcdir)/am/global-rules
