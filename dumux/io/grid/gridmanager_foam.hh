// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Grids
 * \brief Grid manager specialization for FoamGrid
 * \anchor foam_grid_manager
 */
#ifndef DUMUX_IO_GRID_MANAGER_FOAM_HH
#define DUMUX_IO_GRID_MANAGER_FOAM_HH

// FoamGrid specific includes
#if HAVE_DUNE_FOAMGRID
#include <array>
#include <bitset>
#include <limits>
#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/dgffoam.hh>
#include <dune/foamgrid/parallel/distribute.hh>
#endif

#ifndef DUMUX_IO_GRID_MANAGER_BASE_HH
#include <dumux/io/grid/gridmanager_base.hh>
#endif

#include <dumux/common/gridcapabilities.hh>

namespace Dumux {

#if HAVE_DUNE_FOAMGRID

namespace Detail {

/*!
 * \brief Stage a spatial redistribution of a (replicated) FoamGrid over the available ranks, using
 *        the decomposition of a coarse Dune::YaspGrid, when running in parallel.
 *
 * Opt-in via the parameter "<group>.Grid.ParallelPartitioning" (any value enables it). The coarse
 * grid resolution and overlap are controlled by "Grid.ParallelPartitionCells" (default 10) and
 * "Grid.ParallelOverlap" (default 1); the network overlap matches the coarse grid's overlap.
 *
 * This only installs the communicator and *stages* the distribution; the actual rebuild happens in
 * the subsequent GridManager loadBalance() (FoamGrid::loadBalance[(dataHandle)]). Running it there
 * lets the standard DGF/gmsh/vtk parameter migration (Dune::GridPtr / DuMux GridData, keyed by the
 * persistent localIdSet id which the rebuild preserves) remap per-element/-vertex data such as a
 * network radius onto the distributed local elements automatically.
 */
template<class Grid>
void prepareFoamGridDistribution(Grid& grid, const std::string& paramGroup)
{
    if (Dune::MPIHelper::instance().size() <= 1)
        return;
    if (!hasParamInGroup(paramGroup, "Grid.ParallelPartitioning"))
        return; // distribution is opt-in

    static constexpr int dimworld = Grid::dimensionworld;
    using ct = typename Grid::ctype;

    // install the parallel communicator (the grid was read replicated with a self communicator)
    grid.setCommunicator(Dune::MPIHelper::getCommunicator());

    // padded bounding box of the network
    Dune::FieldVector<ct, dimworld> lo(std::numeric_limits<ct>::max());
    Dune::FieldVector<ct, dimworld> hi(std::numeric_limits<ct>::lowest());
    for (const auto& v : vertices(grid.leafGridView()))
    {
        const auto x = v.geometry().center();
        for (int d = 0; d < dimworld; ++d) { lo[d] = std::min(lo[d], x[d]); hi[d] = std::max(hi[d], x[d]); }
    }
    for (int d = 0; d < dimworld; ++d) { const ct m = 0.05*(hi[d]-lo[d]) + 1e-12; lo[d] -= m; hi[d] += m; }

    // coarse host grid that defines the spatial domain decomposition
    const int cells = getParamFromGroup<int>(paramGroup, "Grid.ParallelPartitionCells", 10);
    const int overlap = getParamFromGroup<int>(paramGroup, "Grid.ParallelOverlap", 1);
    using Coarse = Dune::YaspGrid<dimworld, Dune::EquidistantOffsetCoordinates<ct, dimworld>>;
    std::array<int, dimworld> cellsArray; cellsArray.fill(cells);
    Dune::Communication<typename Dune::MPIHelper::MPICommunicator> comm(Dune::MPIHelper::getCommunicator());
    Coarse coarse(lo, hi, cellsArray, std::bitset<dimworld>(0), overlap, comm);

    // compute the partition now (grid still full) and stage the rebuild for loadBalance()
    Dune::FoamGridParallel::prepareSpatialRedistribution(grid, coarse.leafGridView());
}

} // end namespace Detail

/*!
 * \ingroup Grids
 * \brief Provides a grid manager for FoamGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.
 *
 * The following keys are recognized:
 * - File : A dgf/msh/vtp/vtu file to load from, type detection by file extension
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 * - Refinement : the number of global refines to perform
 * - Verbosity : whether the grid construction should output to standard out
 *
 * Additionally if reading from a msh file:
 * - BoundarySegments : whether to insert boundary segments into the grid
 * - DomainMarkers : weather domain marker information should be read from the grid file
 */
template<int dim, int dimworld>
class GridManager<Dune::FoamGrid<dim, dimworld>>
: public GridManagerBase<Dune::FoamGrid<dim, dimworld>>
{
public:
    using Grid = Dune::FoamGrid<dim, dimworld>;
    using ParentType = GridManagerBase<Grid>;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    void init(const std::string& modelParamGroup = "")
    {
        // First, try to create it from file
        if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"), modelParamGroup);
            ParentType::maybeRefineGrid(modelParamGroup);
            Detail::prepareFoamGridDistribution(this->grid(), modelParamGroup);
            ParentType::loadBalance();
            return;
        }

        // Then look for the necessary keys to construct a structured grid from the input file
        else if (hasParamInGroup(modelParamGroup, "Grid.UpperRight"))
        {
            ParentType::template makeStructuredGrid<dim, dimworld>(ParentType::CellType::Simplex, modelParamGroup);
            ParentType::maybeRefineGrid(modelParamGroup);
            Detail::prepareFoamGridDistribution(this->grid(), modelParamGroup);
            ParentType::loadBalance();
        }

        // Didn't find a way to construct the grid
        else
        {
            const auto prefix = modelParamGroup.empty() ? modelParamGroup : modelParamGroup + ".";
            DUNE_THROW(ParameterException, "Please supply one of the parameters "
                                           << prefix + "Grid.UpperRight"
                                           << ", or a grid file in " << prefix + "Grid.File");

        }
    }
};

/*!
 * \ingroup Grids
 * \brief Provides a grid manager for FoamGrids of dim 1
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.

 * The following keys are recognized:
 * - File : A dgf/msh/vtp/vtu file to load from, type detection by file extension
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 * - Refinement : the number of global refines to perform
 * - Verbosity : whether the grid construction should output to standard out
 *
 * Additionally if reading from a msh file:
 * - BoundarySegments : whether to insert boundary segments into the grid
 * - DomainMarkers : weather domain marker information should be read from the grid file
 */
template<int dimworld>
class GridManager<Dune::FoamGrid<1, dimworld>>
: public GridManagerBase<Dune::FoamGrid<1, dimworld>>
{
public:
    using Grid = Dune::FoamGrid<1, dimworld>;
    using ParentType = GridManagerBase<Grid>;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    void init(const std::string& modelParamGroup = "")
    {
        // Try to create it from file
        if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"), modelParamGroup);
            ParentType::maybeRefineGrid(modelParamGroup);
            Detail::prepareFoamGridDistribution(this->grid(), modelParamGroup);
            ParentType::loadBalance();
            return;
        }

        // The required parameters
        using GlobalPosition = Dune::FieldVector<typename Grid::ctype, dimworld>;
        const auto upperRight = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.UpperRight");
        const auto lowerLeft = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.LowerLeft", GlobalPosition(0.0));
        using CellArray = std::array<unsigned int, 1>;
        const auto cells = getParamFromGroup<CellArray>(modelParamGroup, "Grid.Cells", std::array<unsigned int, 1>{{1}});

        // make the grid (structured interval grid in dimworld space)
        Dune::GridFactory<Grid> factory;

        constexpr auto geomType = Dune::GeometryTypes::line;

        // create a step vector
        GlobalPosition step = upperRight;
        step -= lowerLeft, step /= cells[0];

        // create the vertices
        GlobalPosition globalPos = lowerLeft;
        for (unsigned int vIdx = 0; vIdx <= cells[0]; vIdx++, globalPos += step)
            factory.insertVertex(globalPos);

        // create the cells
        for(unsigned int eIdx = 0; eIdx < cells[0]; eIdx++)
            factory.insertElement(geomType, {eIdx, eIdx+1});

        ParentType::gridPtr() = std::shared_ptr<Grid>(factory.createGrid());
        ParentType::maybeRefineGrid(modelParamGroup);
        Detail::prepareFoamGridDistribution(this->grid(), modelParamGroup);
        ParentType::loadBalance();
    }
};

namespace Grid::Capabilities {

// To the best of our knowledge FoamGrid is view thread-safe
// This specialization can be removed after we depend on Dune release 2.9 in which this is guaranteed by FoamGrid itself
template<int dim, int dimworld>
struct MultithreadingSupported<Dune::FoamGrid<dim, dimworld>>
{
    template<class GV>
    static bool eval(const GV&) // default is independent of the grid view
    { return true; }
};

} // end namespace Grid::Capabilities

#endif // HAVE_DUNE_FOAMGRID

} // end namespace Dumux

#endif
