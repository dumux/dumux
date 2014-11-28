check_PROGRAMS = test_box3p3c \
                 test_box3p3cnicolumnxylol \
                 test_box3p3cnikuevette \
                 test_cc3p3c \
                 test_cc3p3cnicolumnxylol \
                 test_cc3p3cnikuevette

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_box3p3c_SOURCES = test_box3p3c.cc
test_box3p3cnicolumnxylol_SOURCES = test_box3p3cnicolumnxylol.cc
test_box3p3cnikuevette_SOURCES = test_box3p3cnikuevette.cc
test_cc3p3c_SOURCES = test_cc3p3c.cc
test_cc3p3cnicolumnxylol_SOURCES = test_cc3p3cnicolumnxylol.cc
test_cc3p3cnikuevette_SOURCES = test_cc3p3cnikuevette.cc

include $(top_srcdir)/am/global-rules
