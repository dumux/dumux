check_PROGRAMS = test_1p

noinst_HEADERS = *.hh
EXTRA_DIST = grids/*.dgf CMakeLists.txt

test_1p_SOURCES = test_1p.cc
test_1p_CXXFLAGS = $(DUNEMPICPPFLAGS)  
test_1p_LDADD = $(DUNEMPILDFLAGS) 

include $(top_srcdir)/am/global-rules
