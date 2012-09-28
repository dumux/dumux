# programs just to build when "make check" is used
check_PROGRAMS = test_adaptive2p2c test_dec2p2c test_multiphysics2p2c

noinst_HEADERS = *.hh

test_adaptive2p2c_SOURCES = test_adaptive2p2c.cc
test_dec2p2c_SOURCES = test_dec2p2c.cc
test_multiphysics2p2c_SOURCES = test_multiphysics2p2c.cc

EXTRA_DIST=*reference.vtu *.input CMakeLists.txt
include $(top_srcdir)/am/global-rules
