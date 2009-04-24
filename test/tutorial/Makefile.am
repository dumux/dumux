# tests where program to build and program to run are equal
NORMALTESTS = tutorial_decoupled tutorial_coupled

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = tutorial_decoupled.cc tutorial_coupled.cc 

tutorial_decoupleddir = $(includedir)/dune/dumux/test/tutorial
tutorial_decoupled_HEADERS = tutorial_soilproperties_coupled.hh tutorial_soilproperties_decoupled.hh tutorialproblem_coupled.hh \
	tutorialproblem_decoupled.hh
tutorial_decoupled_SOURCES = tutorial_decoupled.cc

tutorial_coupled_SOURCES = tutorial_coupled.cc

include $(top_srcdir)/am/global-rules
