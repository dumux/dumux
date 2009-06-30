# setting like in dune-web
BASEDIR=../..
CURDIR=doc/doxygen
WHAT=

doxygendir = $(docdir)/html
doxygen_DATA = $(WHAT)

#we need this in the distribution
EXTRA_DIST = $(WHAT) $(DOXYGENDISTFILES) Doxyfile

all: $(WHAT) $(DOXYGENFILES)

include $(top_srcdir)/am/doxygen
include $(top_srcdir)/am/global-rules

# file to clean only in svn and not in tarball tree
SVNCLEANFILES = html $(WHAT) $(DOXYGENHEADER)

clean-local:
	if test "x$(VPATH)" != "x" && test ! -e $(srcdir)/doxygen-tag; then \
	  if test -e $(top_srcdir)/doc/doxygen/Doxydep; then \
	    rm -rf $(SVNCLEANFILES); \
	  fi; \
	fi
