# tests where program to build and program to run are equal
NORMALTESTS = test_pardiso test_matrixroot

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_pardiso.cc

test_pardiso_SOURCES = test_pardiso.cc

test_matrixroot_SOURCES = test_matrixroot.cc
test_matrixroot_CXXFLAGS=-I/usr/include/gsl 
test_matrixroot_LDADD=-lgsl

include $(top_srcdir)/am/global-rules
