check_PROGRAMS = test_boxco2 \
				 test_ccco2 \
				 test_boxco2ni \
				 test_ccco2ni

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST = CMakeLists.txt $(wildcard *.input) $(wildcard grids/*.dgf) co2values.inc

test_boxco2_SOURCES = test_boxco2.cc
test_ccco2_SOURCES = test_ccco2.cc
test_boxco2ni_SOURCES = test_boxco2ni.cc
test_ccco2ni_SOURCES = test_ccco2ni.cc

include $(top_srcdir)/am/global-rules
