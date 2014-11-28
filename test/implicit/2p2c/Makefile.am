check_PROGRAMS = test_box2p2c \
                 test_box2p2cni \
                 test_cc2p2c \
                 test_cc2p2cni

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_box2p2c_SOURCES = test_box2p2c.cc
test_box2p2cni_SOURCES = test_box2p2cni.cc
test_cc2p2c_SOURCES = test_cc2p2c.cc
test_cc2p2cni_SOURCES = test_cc2p2cni.cc

include $(top_srcdir)/am/global-rules
