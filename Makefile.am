# to the tarball, add all files required to make dunecontrol work
# nicely, the files for the CMake build system, some README files and
# a set of default .opts files for dunecontrol.
EXTRA_DIST = dune.module \
             stamp-vc \
             CMake/Modules/*.cmake \
	     CMakeLists.txt \
	     CHANGELOG CONTRIBUTORS INSTALL.cmake LICENSE \
             debug.opts optim.opts 

SUBDIRS = doc dumux m4 test tutorial

# don't follow the full GNU-standard
# we need automake 1.9 or later
AUTOMAKE_OPTIONS = foreign 1.9

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS= CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

