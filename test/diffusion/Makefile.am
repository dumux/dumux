# tests where program to build and program to run are equal
NORMALTESTS = test_diffusion test_equilibrated

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_diffusion_SOURCES = test_diffusion.cc

test_equilibrated_SOURCES = test_equilibrated.cc

include $(top_srcdir)/am/global-rules
