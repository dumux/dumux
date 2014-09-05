# programs just to build when "make check" is used
check_PROGRAMS = tutorial_decoupled tutorial_coupled

noinst_HEADERS := $(wildcard *.hh)

tutorialdir = $(datadir)/tutorial

tutorial_decoupled_SOURCES = tutorial_decoupled.cc

tutorial_coupled_SOURCES = tutorial_coupled.cc

EXTRA_DIST = \
  CMakeLists.txt $(wildcard *.input) \
  solutions_coupled/$(wildcard *.hh) \
  solutions_coupled/$(wildcard *.diff) \
  solutions_coupled/$(wildcard *.input) \
  solutions_decoupled/$(wildcard *.hh) \
  solutions_decoupled/$(wildcard *.diff) \
  solutions_decoupled/$(wildcard *.input)

include $(top_srcdir)/am/global-rules
