check_PROGRAMS = test_boxmpnc test_ccmpnc test_forchheimer2p test_forchheimer1p test_boxmpnckinetic test_boxmpncthermalnonequil
noinst_HEADERS := $(wildcard *.hh)

EXTRA_DIST:= $(wildcard *.input) CMakeLists.txt $(wildcard grids/*.dgf)
gridsdir=$(datadir)/dumux/grids
grids_DATA=grids/$(wildcard *.dgf)

test_boxmpnc_SOURCES = test_boxmpnc.cc
test_ccmpnc_SOURCES = test_ccmpnc.cc
test_forchheimer2p_SOURCES = test_forchheimer2p.cc
test_forchheimer1p_SOURCES = test_forchheimer1p.cc
test_boxmpnckinetic_SOURCES = test_boxmpnckinetic.cc
test_boxmpncthermalnonequil_SOURCES = test_boxmpncthermalnonequil.cc



include $(top_srcdir)/am/global-rules
