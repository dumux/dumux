# tests where program to build and program to run are equal
NORMALTESTS = tissue_tumor new_tissue_tumor

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

tissue_tumor_SOURCES = tissue_tumor.cc
new_tissue_tumor_SOURCES = new_tissue_tumor.cc
new_tissue_tumor_CXXFLAGS = $(MPI_CPPFLAGS) 
new_tissue_tumor_LDADD = $(MPI_LDFLAGS) 

include $(top_srcdir)/am/global-rules
