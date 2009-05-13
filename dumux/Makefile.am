SUBDIRS = diffusion \
	fractionalflow \
	fvgeometry \
	io \
	material \
	nonlinear \
	onedinndgrid \
	pardiso \
	transport \
	.

# the standard debug streams are put into the libdune
noinst_LTLIBRARIES = libdumux.la
SUBLIBS = onedinndgrid/libonedinndgrid.la

libdumux_la_SOURCES =
libdumux_la_LIBADD = $(SUBLIBS)

dumuxdir = $(includedir)/dune/dumux

include $(top_srcdir)/am/global-rules
