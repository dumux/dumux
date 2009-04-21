# tests where program to build and program to run are equal
NORMALTESTS = new_test_lens test_twophase test_twophase_pnsw test_dmt

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = new_test_lens.cc test_twophase.cc test_twophase_pnsw.cc test_dmt.cc

new_test_lens_SOURCES = new_test_lens.cc

test_twophase_SOURCES = test_twophase.cc

test_twophase_pnsw_SOURCES = test_twophase_pnsw.cc

test_dmt_SOURCES = test_dmt.cc

include $(top_srcdir)/am/global-rules
