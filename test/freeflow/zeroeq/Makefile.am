check_PROGRAMS = test_zeroeq \
                 test_zeroeq_channel

test_zeroeq_SOURCES = test_zeroeq.cc
test_zeroeq_channel_SOURCES = test_zeroeq_channel.cc

EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt
noinst_HEADERS:=$(wildcard *.hh)
include $(top_srcdir)/am/global-rules
