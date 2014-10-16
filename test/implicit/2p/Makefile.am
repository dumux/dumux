check_PROGRAMS = test_box2p test_cc2p

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_box2p_SOURCES = test_box2p.cc
test_cc2p_SOURCES = test_cc2p.cc

include $(top_srcdir)/am/global-rules
