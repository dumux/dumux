# programs just to build when "make check" is used
check_PROGRAMS = tutorial_decoupled tutorial_coupled

noinst_HEADERS := $(wildcard *.hh)

tutorialdir = $(datadir)/tutorial

tutorial_decoupled_SOURCES = tutorial_decoupled.cc

tutorial_coupled_SOURCES = tutorial_coupled.cc

EXTRA_DIST = \
  CMakeLists.txt $(wildcard *.input) \
  $(wildcard solutions_coupled/*.hh) \
  $(wildcard solutions_coupled/*.diff) \
  $(wildcard solutions_coupled/*.input) \
  $(wildcard solutions_decoupled/*.hh) \
  $(wildcard solutions_decoupled/*.diff) \
  $(wildcard solutions_decoupled/*.input)

include $(top_srcdir)/am/global-rules
