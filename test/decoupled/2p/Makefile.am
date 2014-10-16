check_PROGRAMS = test_impes \
                 test_impesadaptive \
                 test_impeswithamg \
                 test_transport \
                 test_mpfa2p \
                 test_3d2p

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_impes_SOURCES = test_impes.cc

test_impesadaptive_SOURCES = test_impesadaptive.cc

test_impeswithamg_SOURCES = test_impeswithamg.cc
test_impeswithamg_CXXFLAGS = -DUSE_AMGBACKEND

test_transport_SOURCES = test_transport.cc

test_mpfa2p_SOURCES = test_mpfa2p.cc

test_3d2p_SOURCES = test_3d2p.cc

include $(top_srcdir)/am/global-rules
