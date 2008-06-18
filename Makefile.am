# $Id: duneproject 4887 2007-04-02 10:06:45Z christi $

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

#SUBDIRS = doc dumux lib test benchmarkproblems
SUBDIRS = doc dumux test benchmarkproblems

#pkgconfigdir = $(libdir)/pkgconfig
#pkgconfig_DATA = dune-mux.pc

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) DOXYGEN="true" CXX="$(CXX)" CC="$(CC)"

#include $(top_srcdir)/am/global-rules
