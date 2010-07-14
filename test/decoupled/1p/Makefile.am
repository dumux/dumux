# tests where program to build and program to run are equal
NORMALTESTS = test_diffusion test_diffusion3d test_mpfa \
			  test_equilibrated test_gradedmesh test_finesubdomainmesh \
			  test_cpgrid test_normals test_dgf 

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_diffusion_SOURCES = test_diffusion.cc

test_mpfa_SOURCES = test_mpfa.cc

test_equilibrated_SOURCES = test_equilibrated.cc

test_gradedmesh_SOURCES = test_gradedmesh.cc

test_finesubdomainmesh_SOURCES = test_finesubdomainmesh.cc

test_dgf_SOURCES = test_dgf.cc

test_diffusion3d_SOURCES = test_diffusion3d.cc
test_diffusion3d_CXXFLAGS = -I../../../dune-cornerpoint
test_diffusion3d_LDADD = ../../../dune-cornerpoint/lib/libdunecornerpoint.la \
                    ../../../dune-cornerpoint/grid/preprocess/libpreprocess.la 

test_cpgrid_SOURCES = test_cpgrid.cc
test_cpgrid_CXXFLAGS = -I../../../dune-cornerpoint
test_cpgrid_LDADD = ../../../dune-cornerpoint/lib/libdunecornerpoint.la \
                    ../../../dune-cornerpoint/grid/preprocess/libpreprocess.la \
			-lblas -llapack

test_normals_SOURCES = test_normals.cc
test_normals_CXXFLAGS = -I../../../dune-cornerpoint
test_normals_LDADD = ../../../dune-cornerpoint/lib/libdunecornerpoint.la \
                    ../../../dune-cornerpoint/grid/preprocess/libpreprocess.la 

include $(top_srcdir)/am/global-rules
