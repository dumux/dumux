# tests where program to build and program to run are equal
NORMALTESTS = dynamic_pc

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = dynamic_pc.cc

dynamic_pc_dir = $(includedir)/dune/dumux/appl/dynamic_pc

dynamic_pc_SOURCES = dynamic_pc.cc

include $(top_srcdir)/am/global-rules
