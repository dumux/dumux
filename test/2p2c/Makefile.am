dist_noinst_DATA = test_2p2c.cc new_test_2p2c.cc new_test_blob.cc

noinst_PROGRAMS = test_2p2c new_test_2p2c new_test_blob

test_2p2c_SOURCES = test_2p2c.cc 
new_test_2p2c_SOURCES = new_test_2p2c.cc
new_test_blob_SOURCES = new_test_blob.cc
new_test_blob_CXXFLAGS = $(MPI_CPPFLAGS) 
new_test_blob_LDADD = $(MPI_LDFLAGS) 

include $(top_srcdir)/am/global-rules
