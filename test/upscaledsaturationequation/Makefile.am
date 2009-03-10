# tests where program to build and program to run are equal
NORMALTESTS = test_upscaledsaturation

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_upscaledsaturation.cc

test_upscaledsaturation_SOURCES = test_upscaledsaturation.cc
test_upscaledsaturation_DEPENDENCIES = $(LOCAL_LIBS)

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
