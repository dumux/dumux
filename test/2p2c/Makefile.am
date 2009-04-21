# tests where program to build and program to run are equal
NORMALTESTS = test_2p2c new_test_2p2c new_test_blob

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_2p2c_SOURCES = test_2p2c.cc 

new_test_2p2c_SOURCES = new_test_2p2c.cc

new_test_blob_SOURCES = new_test_blob.cc
#new_test_blob_CXXFLAGS = $(MPI_CPPFLAGS) 
#new_test_blob_LDADD = $(MPI_LDFLAGS) 

include $(top_srcdir)/am/global-rules
