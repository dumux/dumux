check_PROGRAMS = test_box3p test_cc3p test_box3pniconvection test_cc3pniconvection test_box3pniconduction test_cc3pniconduction

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_box3p_SOURCES = test_box3p.cc
test_cc3p_SOURCES = test_cc3p.cc
test_box3pniconvection_SOURCES = test_box3pniconvection.cc
test_cc3pniconvection_SOURCES = test_cc3pniconvection.cc
test_box3pniconduction_SOURCES = test_box3pniconduction.cc
test_cc3pniconduction_SOURCES = test_box3pniconduction.cc

include $(top_srcdir)/am/global-rules
