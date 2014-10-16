check_PROGRAMS = test_boxco2 test_ccco2

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST = CMakeLists.txt $(wildcard *.input) $(wildcard grids/*.dgf) co2values.inc

test_boxco2_SOURCES = test_boxco2.cc
test_ccco2_SOURCES = test_ccco2.cc

include $(top_srcdir)/am/global-rules
