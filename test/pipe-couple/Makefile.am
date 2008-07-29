# $Id: duneproject 4887 2007-04-02 10:06:45Z christi $

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)

noinst_PROGRAMS = pipe_box_simple_couple

pipe_box_simple_couple_SOURCES = pipe_box_simple_couple.cc

pipe_box_simple_couple_CXXFLAGS = $(MPI_CPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) -I$(top_srcdir)/../dune-mux
pipe_box_simple_couple_LDADD = $(MPI_LDFLAGS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) -L$(top_srcdir)/../dune-mux/dumux -ldumux

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"
