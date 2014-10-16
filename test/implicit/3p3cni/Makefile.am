check_PROGRAMS = test_box3p3cni test_cc3p3cni

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_box3p3cni_SOURCES = test_box3p3cni.cc
test_cc3p3cni_SOURCES = test_cc3p3cni.cc

include $(top_srcdir)/am/global-rules
