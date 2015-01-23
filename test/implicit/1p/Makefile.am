check_PROGRAMS = test_box1p \
                 test_box1p1d3d \
                 test_box1p2d3d \
                 test_box1pniconduction \
                 test_box1pniconvection \
                 test_box1pwithamg \
                 test_cc1p \
                 test_cc1p1d3d \
                 test_cc1p2d3d \
                 test_cc1pniconduction \
                 test_cc1pniconvection \
                 test_cc1pwithamg

noinst_HEADERS := $(wildcard *.hh)
EXTRA_DIST:=$(wildcard *.input) $(wildcard grids/*.dgf) CMakeLists.txt

test_box1p_SOURCES = test_box1p.cc
test_box1p1d3d_SOURCES = test_box1p1d3d.cc
test_box1p2d3d_SOURCES = test_box1p2d3d.cc
test_box1pniconduction_SOURCES = test_box1pniconduction.cc
test_box1pniconvection_SOURCES = test_box1pniconvection.cc
test_box1pwithamg_SOURCES = test_box1pwithamg.cc
test_cc1p_SOURCES = test_cc1p.cc
test_cc1p1d3d_SOURCES = test_cc1p1d3d.cc
test_cc1p2d3d_SOURCES = test_cc1p2d3d.cc
test_cc1pniconduction_SOURCES = test_cc1pniconduction.cc
test_cc1pniconvection_SOURCES = test_cc1pniconvection.cc
test_cc1pwithamg_SOURCES = test_cc1pwithamg.cc

test_cc1p1d3d_LDADD = $(AM_LDFLAGS)  $(ALBERTA3D_LIBS)
test_cc1p1d3d_LDFLAGS= $(AM_LDFLAGS) $(ALBERTA3D_LDFLAGS)
test_cc1p1d3d_CPPFLAGS=$(AM_CPPFLAGS) $(ALBERTA3D_CPPFLAGS)

test_box1p1d3d_LDADD = $(AM_LDFLAGS)  $(ALBERTA3D_LIBS)
test_box1p1d3d_LDFLAGS= $(AM_LDFLAGS) $(ALBERTA3D_LDFLAGS)
test_box1p1d3d_CPPFLAGS=$(AM_CPPFLAGS) $(ALBERTA3D_CPPFLAGS)

test_cc1p2d3d_LDADD = $(AM_LDFLAGS)  $(ALBERTA3D_LIBS)
test_cc1p2d3d_LDFLAGS= $(AM_LDFLAGS) $(ALBERTA3D_LDFLAGS)
test_cc1p2d3d_CPPFLAGS=$(AM_CPPFLAGS) $(ALBERTA3D_CPPFLAGS)

test_box1p2d3d_LDADD = $(AM_LDFLAGS)  $(ALBERTA3D_LIBS)
test_box1p2d3d_LDFLAGS= $(AM_LDFLAGS) $(ALBERTA3D_LDFLAGS)
test_box1p2d3d_CPPFLAGS=$(AM_CPPFLAGS) $(ALBERTA3D_CPPFLAGS)


include $(top_srcdir)/am/global-rules
