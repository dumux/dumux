SUBDIRS = 2p2c 2p2cni 2pni brinkman coupled diffusion fractionalflow functions fvgeometry io material minc nonlinear onedinndgrid operators parabolic pardiso shapefunctions simset stokes timedisc transport twophase .

# the standard debug streams are put into the libdune
noinst_LTLIBRARIES = libdumux.la
SUBLIBS = shapefunctions/libshapefunctions.la onedinndgrid/libonedinndgrid.la

libdumux_la_SOURCES =
libdumux_la_LIBADD = $(SUBLIBS)

dumuxdir = $(includedir)/dune/dumux

include $(top_srcdir)/am/global-rules
