# checks
include(CheckConstexpr)
# additional macros
include(DumuxTestMacros)

find_package(Gnuplot)
set(HAVE_GNUPLOT ${GNUPLOT_FOUND})
