bin_PROGRAMS = test_richards
noinst_HEADERS = *.hh
dist_noinst_DATA = *.cc

EXTRA_DIST=grids/*.dgf
gridsdir=$(datadir)/dumux/grids
grids_DATA=grids/*.dgf

test_richards_SOURCES = test_richards.cc richardsproblem.hh richardssoil.hh

include $(top_srcdir)/am/global-rules
