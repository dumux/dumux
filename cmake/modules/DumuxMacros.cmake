# checks
include(CheckConstexpr)
# additional macros
include(AddInputFileLinks)
include(DumuxTestMacros)

find_package(Gnuplot)
set(HAVE_GNUPLOT ${GNUPLOT_FOUND})
