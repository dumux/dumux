check_PROGRAMS = test_boxrichards test_ccrichards

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_boxrichards_SOURCES = test_boxrichards.cc 
test_ccrichards_SOURCES = test_ccrichards.cc 

include $(top_srcdir)/am/global-rules
