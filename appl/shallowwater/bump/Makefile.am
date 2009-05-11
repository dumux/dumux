# tests where program to build and program to run are equal
NORMALTESTS = 2d_bump_subcritical 2d_bump_transcritical 2d_bump_transcritical_shock

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

2d_bump_subcritical_SOURCES = 2d_bump_subcritical.cc

2d_bump_transcritical_SOURCES = 2d_bump_transcritical.cc

2d_bump_transcritical_shock_SOURCES = 2d_bump_transcritical_shock.cc


include $(top_srcdir)/am/global-rules

