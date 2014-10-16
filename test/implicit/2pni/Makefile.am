check_PROGRAMS = test_box2pni test_cc2pni

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_box2pni_SOURCES = test_box2pni.cc
test_cc2pni_SOURCES = test_cc2pni.cc

include $(top_srcdir)/am/global-rules
