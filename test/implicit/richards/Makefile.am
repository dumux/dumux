check_PROGRAMS = \
	test_boxrichards \
	test_ccrichards \
	test_boxrichardsniconvection \ 
	test_ccrichardsniconvection \ 
	test_boxrichardsniconduction \
	test_ccrichardsniconduction

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_boxrichards_SOURCES = test_boxrichards.cc 
test_ccrichards_SOURCES = test_ccrichards.cc 
test_boxrichardsniconvection_SOURCES = test_boxrichardsniconvection.cc
test_ccrichardsniconvection_SOURCES = test_ccrichardsniconvection.cc
test_boxrichardsniconduction_SOURCES = test_boxrichardsniconduction.cc
test_ccrichardsniconduction_SOURCES = test_ccrichardsniconduction.cc

include $(top_srcdir)/am/global-rules
