dist_noinst_DATA = test_upscaledsaturation.cc test_heterogeneoussoil.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_upscaledsaturation test_heterogeneoussoil

test_upscaledsaturation_SOURCES = test_upscaledsaturation.cc
test_upscaledsaturation_DEPENDENCIES = $(LOCAL_LIBS)
test_upscaledsaturation_CXXFLAGS = $(MPI_CPPFLAGS) 
test_upscaledsaturation_LDADD = $(MPI_LDFLAGS) 

test_heterogeneoussoil_SOURCES = test_heterogeneoussoil.cc
test_heterogeneoussoil_DEPENDENCIES = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
