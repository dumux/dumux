# $Id: duneproject 5842 2010-01-20 18:48:34Z joe $

# we need the module file to be able to build via dunecontrol
EXTRA_DIST = dune.module \
             CMake/Modules/*.cmake CMakeLists.txt \
             CONTRIBUTORS INSTALL.cmake LICENSE \ 
             debug.opts optim.opts 

SUBDIRS = appl doc dumux m4 test tutorial

if BUILD_DOCS
# TODO: set up documentation tree automatically
#SUBDIRS += doc
endif

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

