# tests where program to build and program to run are equal
NORMALTESTS = test_transport test_boxtransport test_chtransport #test_implicitfvtransport

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_transport.cc test_boxtransport.cc 

test_transportdir = $(includedir)/dune/dumux/test/transport
test_transport_HEADERS = simplerproblem.hh
test_transport_SOURCES = test_transport.cc

test_boxtransport_SOURCES = test_boxtransport.cc
test_chtransport_SOURCES = test_chtransport.cc 
test_implicitfvtransport_SOURCES = test_implicitfvtransport.cc

include $(top_srcdir)/am/global-rules
