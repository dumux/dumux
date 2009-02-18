dist_noinst_DATA = test_coupled.cc test_stokesdarcy.cc test_boxstokesdarcy.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_coupled test_stokesdarcy test_boxstokesdarcy

test_coupleddir = $(includedir)/dune/dumux/test/coupled
test_coupled_HEADERS = boxdiffusion.hh boxdiffusionjacobian.hh diffusionparameters.hh 
test_coupled_SOURCES = test_coupled.cc

test_stokesdarcy_SOURCES = test_stokesdarcy.cc

test_boxstokesdarcy_SOURCES = test_boxstokesdarcy.cc

include $(top_srcdir)/am/global-rules
