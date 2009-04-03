# tests where program to build and program to run are equal
NORMALTESTS = injection_ff

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = injection_ff.cc

injection_ffdir = $(includedir)/dune/dumux/appl/phasepressure_ff

injection_ff_SOURCES = injection_ff.cc

include $(top_srcdir)/am/global-rules
