if BUILD_HANDBOOK
TEXSOURCE := $(wildcard *.tex)
MAINDOCS = dumux-handbook.pdf

DOCS = $(HTDOCS) $(MAINDOCS)

# setting like in dune-web
BASEDIR=../..
CURDIR=doc/handbook

handbookdir=$(docdir)/handbook
handbook_DATA = dumux-handbook.pdf

CLEANFILES := $(wildcard *.aux) $(wildcard *.bbl) $(wildcard *.blg) $(wildcard *.log) $(wildcard *.out) $(wildcard *.toc) $(wildcard *.dvi) $(wildcard *.ps) \
  $(wildcard *.4ct) $(wildcard *.4tc) $(wildcard *.css) $(wildcard *.idv) $(wildcard *.idx) $(wildcard *.lg) $(wildcard *.tid) $(wildcard *.tmp) $(wildcard *.tms) $(wildcard *.xref)
CLEANFILES += $(DOCS)

doc: $(TEXSOURCE)
	@PDFLATEX@ dumux-handbook.tex <&-
	@BIBTEX@ dumux-handbook <&-
	@PDFLATEX@ dumux-handbook.tex <&-
	@PDFLATEX@ dumux-handbook.tex <&-

#dist-hook:
#	sed $(srcdir)/Makefile.dist.am -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.am
#	sed $(srcdir)/Makefile.dist.in -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.in

EXTRA_TEXINPUTS=$(top_srcdir)

else
doc:
	echo -n ""
endif

dumux-handbook.pdf:
	$(MAKE) $(AM_MAKEFLAGS) doc

EXTRA_DIST= dumux-handbook.pdf \
		CMakeLists.txt	

include $(top_srcdir)/am/global-rules
