# $Id$

# setting like in dune-web
BASEDIR=../..
CURDIR=doc/doxygen
WHAT=

doxygendir = $(datadir)/doc/dune-grid/doxygen
doxygen_DATA = $(WHAT)

#we need this in the distribution
EXTRA_DIST = $(WHAT) $(DOXYGENDISTFILES) doxy-footer.html doxy-header.wml mainpage Doxydep Doxyfile

all: $(WHAT) $(DOXYGENFILES)

include $(top_srcdir)/am/doxygen
include $(top_srcdir)/am/global-rules
