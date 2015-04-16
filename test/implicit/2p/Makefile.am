check_PROGRAMS = test_box2p \
                 test_box2pni \
                 test_cc2p \
                 test_cc2pni \
                 test_ccadaptive2p

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_box2p_SOURCES = test_box2p.cc
test_box2pni_SOURCES = test_box2pni.cc
test_cc2p_SOURCES = test_cc2p.cc
test_cc2pni_SOURCES = test_cc2pni.cc
test_ccadaptive2p_SOURCES = test_ccadaptive2p.cc

include $(top_srcdir)/am/global-rules
