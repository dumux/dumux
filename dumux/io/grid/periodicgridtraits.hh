// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Grid properties related to periodicity
 */

#ifndef DUMUX_IO_GRID_PERIODIC_GRID_TRAITS_HH
#define DUMUX_IO_GRID_PERIODIC_GRID_TRAITS_HH

#include <type_traits>
#include <dune/common/exceptions.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/grid/common/exceptions.hh>

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#endif

#if HAVE_DUNE_SUBGRID
#include <dune/subgrid/subgrid.hh>
#endif

namespace Dumux {

template<typename Grid>
struct PeriodicGridTraits
{
    struct SupportsPeriodicity : public std::false_type {};

    PeriodicGridTraits(const Grid& grid) {};

    bool isPeriodic (const typename Grid::LeafIntersection& intersection) const
    {
        return false;
    }
};

#if HAVE_DUNE_SPGRID
template<class ct, int dim, template< int > class Ref, class Comm>
struct PeriodicGridTraits<Dune::SPGrid<ct, dim, Ref, Comm>>
{
private:
    using Grid = Dune::SPGrid<ct, dim, Ref, Comm>;
public:
    struct SupportsPeriodicity : public std::true_type {};

    PeriodicGridTraits(const Grid& grid) {};

    bool isPeriodic (const typename Grid::LeafIntersection& intersection) const
    {
        return intersection.neighbor() && intersection.boundary();
    }
};
#endif //HAVE_DUNE_SPGRID

#if HAVE_DUNE_SUBGRID
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

    bool isPeriodic (const typename Grid::LeafIntersection& intersection) const
    {
        const auto& hostElement = subGrid_.template getHostEntity<0>(intersection.inside());
        for (const auto& hostIntersection : intersections(subGrid_.getHostGrid().leafGridView(), hostElement))
        {
            if (hostIntersection.indexInInside() == intersection.indexInInside())
            {
                const bool periodicInHostGrid = hostTraits_.isPeriodic(hostIntersection);
                return periodicInHostGrid && subGrid_.template contains<0>(hostIntersection.outside());
            }
        }
        return false;
    }

    void verifyConformingPeriodicBoundary() const
    {
        for (const auto& element : elements(subGrid_.leafGridView()))
        {
            for (const auto& intersection : intersections(subGrid_.leafGridView(), element))
            {
                const auto& hostElement = subGrid_.template getHostEntity<0>(intersection.inside());
                for (const auto& hostIntersection : intersections(subGrid_.getHostGrid().leafGridView(), hostElement))
                {
                    if (hostIntersection.indexInInside() == intersection.indexInInside())
                    {
                        const bool periodicInHostGrid = hostTraits_.isPeriodic(hostIntersection);
                        if (periodicInHostGrid && !subGrid_.template contains<0>(hostIntersection.outside()))
                            DUNE_THROW(Dune::GridError, "Periodic boundary in host grid but outside"
                                    << " element not included in subgrid. If this is intentional,"
                                    << " take additional care with boundary conditions and remove"
                                    << " verification call.");
                        break;
                    }
                }
            }
        }
    }
};
#endif //HAVE_DUNE_SUBGRID

template<class T>
class SupportsPeriodicity
{
    template<class G>
    using SP = typename G::SupportsPeriodicity;
public:
    using type = typename Dune::Std::detected_or<std::false_type, SP, T>::type;
};

template<class T>
static constexpr bool supportsPeriodicity()
{ return typename SupportsPeriodicity<T>::type(); }

} // end namespace Dumux

#endif
