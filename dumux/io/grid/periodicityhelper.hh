
#ifndef DUMUX_DISCRETIZATION_PERIODICITY_HELPER_HH
#define DUMUX_DISCRETIZATION_PERIODICITY_HELPER_HH

#include <dune/grid/spgrid.hh>

namespace Dumux {

template<typename Grid>
struct PeriodicityHelper
{
    template<typename Intersection>
    static bool isPeriodic (const Intersection& intersection)
    {
        return false;
    }
};

template<class ct, int dim, template< int > class Ref, class Comm>
struct PeriodicityHelper<Dune::SPGrid<ct, dim, Ref, Comm>>
{
    template<typename Intersection>
    static bool isPeriodic (const Intersection& intersection)
    {
        return intersection.neighbor() && intersection.boundary();
    }
};

} // end namespace Dumux

#endif
