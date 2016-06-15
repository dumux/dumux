# additional macros
include(AddInputFileLinks)
include(DumuxDoxygen)
include(DumuxTestMacros)

find_package(Gnuplot)
set(HAVE_GNUPLOT ${GNUPLOT_FOUND})

find_package(Valgrind)
find_package(GLPK)
find_package(NLOPT)
