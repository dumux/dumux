SUBDIRS = diffusion \
	fractionalflow \
	fvgeometry \
	io \
	material \
	nonlinear \
	onedinndgrid \
	pardiso \
	shapefunctions \
	transport \
	.

# the standard debug streams are put into the libdune
noinst_LTLIBRARIES = libdumux.la
SUBLIBS = shapefunctions/libshapefunctions.la onedinndgrid/libonedinndgrid.la

libdumux_la_SOURCES =
libdumux_la_LIBADD = $(SUBLIBS)

dumuxdir = $(includedir)/dune/dumux

include $(top_srcdir)/am/global-rules
