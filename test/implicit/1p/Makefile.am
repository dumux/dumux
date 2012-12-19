check_PROGRAMS = test_box1p test_cc1p

noinst_HEADERS = *.hh
EXTRA_DIST=*reference.vtu *.input grids/*.dgf CMakeLists.txt

test_box1p_SOURCES = test_box1p.cc
test_box1p_CXXFLAGS = $(DUNEMPICPPFLAGS)  
test_box1p_LDADD = $(DUNEMPILDFLAGS) 

test_cc1p_SOURCES = test_cc1p.cc
test_cc1p_CXXFLAGS = $(DUNEMPICPPFLAGS)  
test_cc1p_LDADD = $(DUNEMPILDFLAGS) 

include $(top_srcdir)/am/global-rules
