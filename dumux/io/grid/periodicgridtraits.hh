
#ifndef DUMUX_IO_GRID_PERIODIC_GRID_TRAITS_HH
#define DUMUX_IO_GRID_PERIODIC_GRID_TRAITS_HH

// forward declare
namespace Dune {

template<class ct, int dim, template< int > class Ref, class Comm>
class SPGrid;

template<int dim, class HostGrid, bool MapIndexStorage>
class SubGrid;
} // end namespace Dune

namespace Dumux {

template<typename Grid>
struct PeriodicGridTraits
{
    struct SupportsPeriodicity : public std::false_type {};

    PeriodicGridTraits() {};

    PeriodicGridTraits(const Grid& grid) {};

    template<typename Intersection>
    bool isPeriodic (const Intersection& intersection) const
    {
        return false;
    }
};

template<class ct, int dim, template< int > class Ref, class Comm>
struct PeriodicGridTraits<Dune::SPGrid<ct, dim, Ref, Comm>>
{
private:
    using Grid = Dune::SPGrid<ct, dim, Ref, Comm>;
public:
    struct SupportsPeriodicity : public std::true_type {};

    PeriodicGridTraits() {};

    PeriodicGridTraits(const Grid& grid) {};

    template<typename Intersection>
    bool isPeriodic (const Intersection& intersection) const
    {
        return intersection.neighbor() && intersection.boundary();
    }
};

// SubGrid does not preserve intersection.boundary() at periodic boundaries of host grid
template<int dim, typename HostGrid, bool MapIndexStorage>
struct PeriodicGridTraits<Dune::SubGrid<dim, HostGrid, MapIndexStorage>>
{
private:
    using Grid = Dune::SubGrid<dim, HostGrid, MapIndexStorage>;

    const Grid& subGrid_;
    const PeriodicGridTraits<HostGrid> hostTraits_;

public:
    struct SupportsPeriodicity : public PeriodicGridTraits<HostGrid>::SupportsPeriodicity {};

    PeriodicGridTraits(const Grid& subGrid)
        : subGrid_(subGrid), hostTraits_(subGrid_.getHostGrid()) {};

    template<typename Intersection>
    bool isPeriodic (const Intersection& intersection) const
    {
        const auto& hostElement = subGrid_.template getHostEntity<0>(intersection.inside());
        for (const auto& hostIntersection : intersections(subGrid_.getHostGrid().leafGridView(), hostElement))
        {
            if (hostIntersection.indexInInside() == intersection.indexInInside())
            {
                const bool periodicInHostGrid = hostTraits_.isPeriodic(hostIntersection);
                if (periodicInHostGrid && !subGrid_.template contains<0>(hostIntersection.outside()))
                    DUNE_THROW(Dune::GridError, "Periodic boundary in host grid but outside element not included in subgrid");
                return periodicInHostGrid;
            }
        }
        return false;
    }
};

} // end namespace Dumux

#endif
