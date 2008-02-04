dist_noinst_DATA = test_fractionalflow.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_fractionalflow test_2p2c_decoupled
test_fractionalflow_SOURCES = test_fractionalflow.cc
test_fractionalflow_DEPENDENCIES = $(LOCAL_LIBS)

test_2p2c_decoupled_SOURCES = test_2p2c_decoupled.cc
test_2p2c_decoupled_DEPENDENCIES = $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
