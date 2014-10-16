check_PROGRAMS = test_box1p test_cc1p test_box1pwithamg test_cc1pwithamg

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_box1p_SOURCES = test_box1p.cc
test_cc1p_SOURCES = test_cc1p.cc
test_box1pwithamg_SOURCES = test_box1pwithamg.cc
test_cc1pwithamg_SOURCES = test_cc1pwithamg.cc

include $(top_srcdir)/am/global-rules
