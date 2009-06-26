# $Id$

# setting like in dune-web
BASEDIR=../..
CURDIR=doc/doxygen

if BUILD_DOCS
htmldir = $(datadir)/doc/dumux/html
WHAT=html/*
html_DATA = $(WHAT)
else
WHAT=
endif

#we need this in the distribution
EXTRA_DIST = $(WHAT) $(DOXYGENDISTFILES) doxy-footer.html doxy-header.wml mainpage Doxydep Doxyfile

all: $(WHAT) $(DOXYGENFILES)

include $(top_srcdir)/am/doxygen
include $(top_srcdir)/am/global-rules
