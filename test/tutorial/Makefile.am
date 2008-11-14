dist_noinst_DATA = tutorial_decoupled.cc tutorial_coupled.cc 

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = tutorial_decoupled tutorial_coupled

tutorial_decoupleddir = $(includedir)/dune/dumux/test/tutorial
tutorial_decoupled_HEADERS = tutorial_soilproperties_coupled.hh tutorial_soilproperties_decoupled.hh tutorialproblem_coupled.hh \
	tutorialproblem_decoupled.hh
tutorial_decoupled_SOURCES = tutorial_decoupled.cc

tutorial_coupled_SOURCES = tutorial_coupled.cc

include $(top_srcdir)/am/global-rules
