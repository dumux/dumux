# to the tarball, add all files required to make dunecontrol work
# nicely, the files for the CMake build system, some README files and
# a set of default .opts files for dunecontrol.
EXTRA_DIST = \
  bin/fuzzycomparevtu.py \
  bin/runTest.sh  \
  cmake/modules/CMakeLists.txt \
  cmake/modules/$(wildcard *.cmake) \
  CMakeLists.txt \
  config.h.cmake \
  INSTALL.cmake \
  dune.module \
  CHANGELOG COPYING README \
  debug.opts optim.opts 

SUBDIRS = doc dumux m4 patches test tutorial

# don't follow the full GNU-standard
# we need automake 1.9 or later
AUTOMAKE_OPTIONS = foreign 1.9

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS= CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

