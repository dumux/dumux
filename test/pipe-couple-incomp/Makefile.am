dist_noinst_DATA = pipe_box_incomp_couple.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = pipe_box_incomp_couple

pipe_box_incomp_coupledir = $(includedir)/dune/dumux/test/pipe-couple-incomp
pipe_box_incomp_couple_HEADERS = dune_multiscale_headers.hh pipephysic.hh properties.hh vtkout.hh pipeflow_incomp_couple.hh boxdiffusion_couple.hh boxdiffusionjacobian_couple.hh diffusionparameters_couple.hh 
pipe_box_incomp_couple_SOURCES = pipe_box_incomp_couple.cc
pipe_box_incomp_couple_DEPENDENCIES = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
