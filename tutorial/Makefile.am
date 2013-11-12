# programs just to build when "make check" is used
check_PROGRAMS = tutorial_decoupled tutorial_coupled

noinst_HEADERS = *.hh

tutorialdir = $(datadir)/tutorial

tutorial_decoupled_SOURCES = tutorial_decoupled.cc

tutorial_coupled_SOURCES = tutorial_coupled.cc

EXTRA_DIST = \
  CMakeLists.txt *.input \
  solutions_coupled/*.hh \
  solutions_coupled/*.diff \
  solutions_coupled/*.input \
  solutions_decoupled/*.hh \
  solutions_decoupled/*.diff \
  solutions_decoupled/*.input

include $(top_srcdir)/am/global-rules
