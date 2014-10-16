# setting like in dune-web
BASEDIR=../..
CURDIR=doc/doxygen

# add some more stuff to install and tarball
DOCFILES = $(PAGES)
EXTRA_DIST = $(PAGES) mainpage modules $(wildcard images/*.png) extradoc/parameterlist extradoc/parameters.html CMakeLists.txt

include $(top_srcdir)/am/doxygen
include $(top_srcdir)/am/global-rules

