SUBDIRS = diffusion fractionalflow functions io material nonlinear onedinndgrid operators parabolic shapefunctions timedisc transport twophase .

# the standard debug streams are put into the libdune
noinst_LTLIBRARIES = libdumux.la
SUBLIBS = shapefunctions/libshapefunctions.la onedinndgrid/libonedinndgrid.la

libdumux_la_SOURCES =
libdumux_la_LIBADD = $(SUBLIBS)

dumuxdir = $(includedir)/dune/dumux

include $(top_srcdir)/am/global-rules
