
#ifndef DUMUX_IO_GRID_PERIODICITY_HELPER_HH
#define DUMUX_IO_GRID_PERIODICITY_HELPER_HH

#include <dune/grid/spgrid.hh>

namespace Dumux {

template<typename Grid>
struct PeriodicityHelper
{
    PeriodicityHelper<Grid>() {};

    PeriodicityHelper<Grid>(const Grid& grid) {};

    template<typename Intersection>
    bool isPeriodic (const Intersection& intersection) const
    {
        return false;
    }
};

template<class ct, int dim, template< int > class Ref, class Comm>
struct PeriodicityHelper<Dune::SPGrid<ct, dim, Ref, Comm>>
{
private:
    using Grid = Dune::SPGrid<ct, dim, Ref, Comm>;
public:
    PeriodicityHelper<Grid>() {};

    PeriodicityHelper<Grid>(const Grid& grid) {};

    template<typename Intersection>
    bool isPeriodic (const Intersection& intersection) const
    {
        return intersection.neighbor() && intersection.boundary();
    }
};

} // end namespace Dumux

#endif
