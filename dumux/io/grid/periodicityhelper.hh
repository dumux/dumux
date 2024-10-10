
#ifndef DUMUX_DISCRETIZATION_PERIODICITY_HELPER_HH
#define DUMUX_DISCRETIZATION_PERIODICITY_HELPER_HH

#include <dune/common/float_cmp.hh>

#include <dune/grid/spgrid.hh>
#include <dune/grid/yaspgrid.hh>

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

template<class ct, int dim>
struct PeriodicityHelper<Dune::YaspGrid<dim, ct>>
{
    template<typename Intersection>
    static bool isPeriodic (const Intersection& intersection)
    {
        if (!intersection.neighbor())
            return false;
        const auto distance = intersection.outside().geometry().center() - intersection.inside().geometry().center();
        const auto& intersectionUnitOuterNormal = intersection.centerUnitOuterNormal();
        const auto reverseNormal = distance/distance.two_norm();
        if (Dune::FloatCmp::eq(intersectionUnitOuterNormal*reverseNormal, -1.0, 1e-7))
            return true;
        return false;
    }
};

template<typename HostIntersection, typename SubIntersection>
bool isEquivalent(const HostIntersection& host, const SubIntersection& sub)
{
    const auto& hCenter = host.geometry().center();
    const auto& sCenter = sub.geometry().center();
    return ((hCenter-sCenter).two_norm2() < 1e-12);
            //&& std::abs(host.unitOuterNormal()*sub.unitOuterNormal()-1.0) < 1e-6);
}

template<int dim, typename HostGrid>
struct PeriodicityHelper<Dune::SubGrid<dim, HostGrid>>
{
    template<typename Intersection, typename SubGrid>
    static bool isPeriodic (const Intersection& intersection, const SubGrid& subGrid)
    {
        const auto& hostElement = subGrid.getHostEntity<0>(subGrid.leafGridView(), intersection.inside());
        for (const auto& hostIntersection : intersections(subGrid.getHostGrid().leafGridView(), hostElement))
        {
            if (isEquivalent(hostIntersection, intersection))
            {
                if (!subGrid.contains(hostIntersection.outside()))
                    return false;
                return PeriodicityHelper<HostGrid>::isPeriodic(hostIntersection);
            }
        }
    }
}

} // end namespace Dumux

#endif
