# tests where program to build and program to run are equal
NORMALTESTS = injection_ff test_decoupled_phasepressure

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = injection_ff.cc test_decoupled_phasepressure.cc

injection_ffdir = $(includedir)/dune/dumux/test/phasepressure_ff

injection_ff_SOURCES = injection_ff.cc

test_decoupled_phasepressure_SOURCES = test_decoupled_phasepressure.cc

include $(top_srcdir)/am/global-rules
