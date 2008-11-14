dist_noinst_DATA = test_box3d.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_box3d

test_box3ddir = $(includedir)/dune/dumux/test/box3d
test_box3d_HEADERS = boxdiffusion.hh boxdiffusionjacobian.hh diffusionparameters.hh
test_box3d_SOURCES = test_box3d.cc

include $(top_srcdir)/am/global-rules
