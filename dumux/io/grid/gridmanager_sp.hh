// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Grid manager specialization for SPGrid
 */
#ifndef DUMUX_IO_GRID_MANAGER_SP_HH
#define DUMUX_IO_GRID_MANAGER_SP_HH

// SPGrid specific includes
#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>
#endif

#ifndef DUMUX_IO_GRID_MANAGER_BASE_HH
#include <dumux/io/grid/gridmanager_base.hh>
#endif

namespace Dumux {

#if HAVE_DUNE_SPGRID

/*!
 * \ingroup InputOutput
 * \brief Provides a grid manager for SPGrid
 *
 * The following keys are recognized:
 * - File : A DGF or gmsh file to load from, type detection by file extension
 *
 */
template<class ct, int dim, template< int > class Ref, class Comm>
class GridManager<Dune::SPGrid<ct, dim, Ref, Comm>>
: public GridManagerBase<Dune::SPGrid<ct, dim, Ref, Comm>>
{
public:
    using Grid = Dune::SPGrid<ct, dim, Ref, Comm>;
    using ParentType = GridManagerBase<Grid>;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    void init(const std::string& paramGroup = "")
    {
        const auto overlap = getParamFromGroup<int>(paramGroup, "Grid.Overlap", 1);
        if (overlap == 0)
            DUNE_THROW(Dune::NotImplemented, "dune-spgrid does currently not support zero overlap!");

        // try to create it from file
        if (hasParamInGroup(paramGroup, "Grid.File"))
        {
            ParentType::makeGridFromDgfFile(getParamFromGroup<std::string>(paramGroup, "Grid.File"));
            ParentType::maybeRefineGrid(paramGroup);
            ParentType::loadBalance();
            return;
        }
        // Didn't find a way to construct the grid
        else if (hasParamInGroup(paramGroup, "Grid.UpperRight"))
        {
            using GlobalPosition = Dune::FieldVector<ct, dim>;
            const auto lowerLeft = getParamFromGroup<GlobalPosition>(paramGroup, "Grid.LowerLeft", GlobalPosition(0.0));
            const auto upperRight = getParamFromGroup<GlobalPosition>(paramGroup, "Grid.UpperRight");

            using IntArray = std::array<int, dim>;
            IntArray cells; cells.fill(1);
            cells = getParamFromGroup<IntArray>(paramGroup, "Grid.Cells", cells);

            const auto periodic = getParamFromGroup<std::bitset<dim>>(paramGroup, "Grid.Periodic", std::bitset<dim>{});
            init(lowerLeft, upperRight, cells, paramGroup, overlap, periodic);
        }
        else
        {
            const auto prefix = paramGroup.empty() ? paramGroup : paramGroup + ".";
            DUNE_THROW(ParameterException, "Please supply a grid file in " << prefix << "Grid.File or " << prefix << "Grid.UpperRight/Cells.");
        }
    }

    void init(const Dune::FieldVector<ct, dim>& lowerLeft,
              const Dune::FieldVector<ct, dim>& upperRight,
              const std::array<int, dim>& cells,
              const std::string& paramGroup = "",
              const int overlap = 1,
              const std::bitset<dim> periodic = std::bitset<dim>{})
    {
        if (overlap == 0)
            DUNE_THROW(Dune::NotImplemented, "dune-spgrid does currently not support zero overlap!");
        using IntArray = std::array<int, dim>;
        IntArray spOverlap; spOverlap.fill(overlap);
        using Domain = typename Grid::Domain;
        std::vector< typename Domain::Cube > cubes;
        cubes.push_back( typename Domain::Cube( lowerLeft, upperRight ) );
        Domain domain( cubes, typename Domain::Topology( static_cast<unsigned int>(periodic.to_ulong()) ) );
        ParentType::gridPtr() = std::make_shared<Grid>( domain, cells, spOverlap );
        ParentType::maybeRefineGrid(paramGroup);
        ParentType::loadBalance();
    }
};

#endif // HAVE_DUNE_SPGRID

} // end namespace Dumux

#endif
