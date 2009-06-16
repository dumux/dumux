SUBDIRS = auxiliary \
        boxmodels \
        diffusion \
	fractionalflow \
        functions \
	fvgeometry \
	functions \
	io \
	material \
	nonlinear \
	onedinndgrid \
	operators \
	pardiso \
	shapefunctions \
        timedisc \
	transport \
	.

# the standard debug streams are put into the libdune
noinst_LTLIBRARIES = libdumux.la
SUBLIBS = shapefunctions/libshapefunctions.la onedinndgrid/libonedinndgrid.la

libdumux_la_SOURCES =
libdumux_la_LIBADD = $(SUBLIBS)

dumuxdir = $(includedir)/dumux
dumux_HEADERS = exceptions.hh 

include $(top_srcdir)/am/global-rules
