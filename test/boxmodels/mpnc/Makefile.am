check_PROGRAMS = test_mpnc test_forchheimer2p test_forchheimer1p
noinst_HEADERS = *.hh

EXTRA_DIST=*reference.vtu *.input CMakeLists.txt grids/*.dgf
gridsdir=$(datadir)/dumux/grids
grids_DATA=grids/*.dgf

test_mpnc_SOURCES = test_mpnc.cc
test_forchheimer2p_SOURCES = test_forchheimer2p.cc
test_forchheimer1p_SOURCES = test_forchheimer1p.cc



include $(top_srcdir)/am/global-rules
