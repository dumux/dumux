check_PROGRAMS = test_box1p2c \
                 test_box1p2cniconduction \
                 test_box1p2cniconvection \
                 test_cc1p2c \
                 test_cc1p2cniconduction \
                 test_cc1p2cniconvection

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_box1p2c_SOURCES = test_box1p2c.cc
test_box1p2cniconduction_SOURCES = test_box1p2cniconduction.cc
test_box1p2cniconvection_SOURCES = test_box1p2cniconvection.cc
test_cc1p2c_SOURCES = test_cc1p2c.cc
test_cc1p2cniconduction_SOURCES = test_cc1p2cniconduction.cc
test_cc1p2cniconvection_SOURCES = test_cc1p2cniconvection.cc

include $(top_srcdir)/am/global-rules
