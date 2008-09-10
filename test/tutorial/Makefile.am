dist_noinst_DATA = tutorial_decoupled.cc tutorial_coupled.cc 

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = tutorial_decoupled tutorial_coupled

tutorial_decoupled_SOURCES = tutorial_decoupled.cc
tutorial_decoupled_DEPENDENCIES = $(LOCAL_LIBS)

tutorial_coupled_SOURCES = tutorial_decoupled.cc
tutorial_coupled_DEPENDENCIES = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
