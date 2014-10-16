check_PROGRAMS = test_boxco2ni test_ccco2ni

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST = CMakeLists.txt $(wildcard *.input) $(wildcard grids/*.dgf) co2values.inc

test_boxco2ni_SOURCES = test_boxco2ni.cc
test_ccco2ni_SOURCES = test_ccco2ni.cc

include $(top_srcdir)/am/global-rules

