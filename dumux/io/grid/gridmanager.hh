// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup InputOutput
 * \brief Provides a grid manager for all supported grid managers with
 *        input file interfaces. Manages data via the grid data member.
 *
 * \todo add Petrel grids with dune-cornerpoint
 */
#ifndef DUMUX_IO_GRID_MANAGER_HH
#define DUMUX_IO_GRID_MANAGER_HH

#include <array>
#include <bitset>
#include <memory>
#include <sstream>

#include <dune/common/exceptions.hh>
#include <dune/common/classname.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

// YaspGrid specific includes
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

 // OneDGrid specific includes
#include <dune/grid/onedgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfoned.hh>

// UGGrid specific includes
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif

// ALUGrid specific includes
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#endif

// FoamGrid specific includes
#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/dgffoam.hh>
#endif

// SPGrid specific includes
#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>
#endif

#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>

#include "griddata.hh"

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief The grid manager base interface (public) and methods common
 *        to most grid manager specializations (protected).
 */
template <class GridType>
class GridManagerBase
{
public:
    using Grid = GridType;
    using GridData = Dumux::GridData<Grid>;

    /*!
     * \brief Make the grid. Implement this method in the specialization of this class for a grid type.
     */
    void init(const std::string& modelParamGroup = "")
    {
        DUNE_THROW(Dune::NotImplemented,
            "The GridManager for grid type " << Dune::className<Grid>() << " is not implemented! Consider providing your own GridManager.");
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    Grid& grid()
    {
        if(enableDgfGridPointer_)
            return *dgfGridPtr();
        else
            return *gridPtr();
    }

    /*!
     * \brief Call loadBalance() function of the grid.
     */
    void loadBalance()
    {
        if (Dune::MPIHelper::getCollectiveCommunication().size() > 1)
        {
            // if we may have dgf parameters use load balancing of the dgf pointer
            if(enableDgfGridPointer_)
            {
                dgfGridPtr().loadBalance();
                // update the grid data
                gridData_ = std::make_shared<GridData>(dgfGridPtr());
            }

            // if we have gmsh parameters we have to manually load balance the data
            else if (enableGmshDomainMarkers_)
            {
                // element and face markers are communicated during load balance
                auto dh = gridData_->createGmshDataHandle();
                gridPtr()->loadBalance(dh.interface());
                gridPtr()->communicate(dh.interface(), Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
            }
            else
                gridPtr()->loadBalance();
        }
    }

    std::shared_ptr<GridData> getGridData() const
    {
        if (!enableDgfGridPointer_ && !enableGmshDomainMarkers_)
            DUNE_THROW(Dune::IOError, "No grid data available");

        return gridData_;
    }


protected:

    /*!
     * \brief Returns a reference to the grid pointer (std::shared_ptr<Grid>)
     */
    std::shared_ptr<Grid>& gridPtr()
    {
        if(!enableDgfGridPointer_)
            return gridPtr_;
        else
            DUNE_THROW(Dune::InvalidStateException, "You are using DGF. To get the grid pointer use method dgfGridPtr()!");
    }

    /*!
     * \brief Returns a reference to the DGF grid pointer (Dune::GridPtr<Grid>).
     */
    Dune::GridPtr<Grid>& dgfGridPtr()
    {
        if(enableDgfGridPointer_)
            return dgfGridPtr_;
        else
            DUNE_THROW(Dune::InvalidStateException, "The DGF grid pointer is only available if the grid was constructed with a DGF file!");
    }

    /*!
     * \brief Returns the filename extension of a given filename
     */
    std::string getFileExtension(const std::string& fileName) const
    {
        std::size_t i = fileName.rfind('.', fileName.length());
        if (i != std::string::npos)
        {
            return(fileName.substr(i+1, fileName.length() - i));
        }
        else
        {
            DUNE_THROW(Dune::IOError, "Please provide and extension for your grid file ('"<< fileName << "')!");
        }
        return "";
    }

    /*!
     * \brief Makes a grid from a file. We currently support *.dgf (Dune Grid Format) and *.msh (Gmsh mesh format).
     */
    void makeGridFromFile(const std::string& fileName,
                          const std::string& modelParamGroup)
    {
        // We found a file in the input file...does it have a supported extension?
        const std::string extension = getFileExtension(fileName);
        if (extension != "dgf" && extension != "msh")
            DUNE_THROW(Dune::IOError, "Grid type " << Dune::className<Grid>() << " only supports DGF (*.dgf) and Gmsh (*.msh) grid files but the specified filename has extension: *."<< extension);

        // make the grid
        if (extension == "dgf")
        {
            enableDgfGridPointer_ = true;
            dgfGridPtr() = Dune::GridPtr<Grid>(fileName.c_str(), Dune::MPIHelper::getCommunicator());
            gridData_ = std::make_shared<GridData>(dgfGridPtr_);
        }
        else if (extension == "msh")
        {
            // get some optional parameters
            const bool verbose = getParamFromGroup<bool>(modelParamGroup, "Grid.Verbosity", false);
            const bool boundarySegments = getParamFromGroup<bool>(modelParamGroup, "Grid.BoundarySegments", false);
            const bool domainMarkers = getParamFromGroup<bool>(modelParamGroup, "Grid.DomainMarkers", false);

            if (domainMarkers)
                enableGmshDomainMarkers_ = true;

            // as default read it on all processes in parallel
            if(domainMarkers)
            {
                std::vector<int> boundaryMarkers, elementMarkers;
                auto gridFactory = std::make_unique<Dune::GridFactory<Grid>>();
                Dune::GmshReader<Grid>::read(*gridFactory, fileName, boundaryMarkers, elementMarkers, verbose, boundarySegments);
                gridPtr() = std::shared_ptr<Grid>(gridFactory->createGrid());
                gridData_ = std::make_shared<GridData>(gridPtr_, std::move(gridFactory), std::move(elementMarkers), std::move(boundaryMarkers));
            }
            else
            {
                auto gridFactory = std::make_unique<Dune::GridFactory<Grid>>();
                Dune::GmshReader<Grid>::read(*gridFactory, fileName, verbose, boundarySegments);
                gridPtr() = std::shared_ptr<Grid>(gridFactory->createGrid());
            }
        }
    }

    /*!
     * \brief Makes a grid from a DGF file. This is used by grid managers that only support DGF.
     */
    void makeGridFromDgfFile(const std::string& fileName)
    {
        // We found a file in the input file...does it have a supported extension?
        const std::string extension = getFileExtension(fileName);
        if(extension != "dgf")
            DUNE_THROW(Dune::IOError, "Grid type " << Dune::className<Grid>() << " only supports DGF (*.dgf) but the specified filename has extension: *."<< extension);

        enableDgfGridPointer_ = true;
        dgfGridPtr() = Dune::GridPtr<Grid>(fileName.c_str(), Dune::MPIHelper::getCommunicator());
        gridData_ = std::make_shared<GridData>(dgfGridPtr_);
    }

    /*!
     * \brief The cell types for structured grids
     */
    enum CellType {Simplex, Cube};

    /*!
     * \brief Makes a structured cube grid using the structured grid factory
     */
    template <int dim, int dimworld>
    void makeStructuredGrid(CellType cellType,
                            const std::string& modelParamGroup)
    {
        using GlobalPosition = Dune::FieldVector<typename Grid::ctype, dimworld>;
        const auto upperRight = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.UpperRight");
        const auto lowerLeft = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.LowerLeft", GlobalPosition(0.0));

        using CellArray = std::array<unsigned int, dim>;
        CellArray cells; cells.fill(1);
        cells = getParamFromGroup<CellArray>(modelParamGroup, "Grid.Cells", cells);

        // make the grid
        if (cellType == CellType::Cube)
        {
            gridPtr() = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, cells);
        }
        else if (cellType == CellType::Simplex)
        {
            gridPtr() = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lowerLeft, upperRight, cells);
        }
        else
        {
            DUNE_THROW(Dune::GridError, "Unknown cell type for making structured grid! Choose Cube or Simplex.");
        }
    }

    /*!
     * \brief Refines a grid after construction if GridParameterGroup.Refinement is set in the input file
     */
    void maybeRefineGrid(const std::string& modelParamGroup)
    {
        if (hasParamInGroup(modelParamGroup, "Grid.Refinement"))
            grid().globalRefine(getParamFromGroup<int>(modelParamGroup, "Grid.Refinement"));
    }

    /*!
    * \brief A state variable if the DGF Dune::GridPtr has been enabled.
    *        It is always enabled if a DGF grid file was used to create the grid.
    */
    bool enableDgfGridPointer_ = false;

    /*!
    * \brief A state variable if domain markers have been read from a Gmsh file.
    */
    bool enableGmshDomainMarkers_ = false;

    std::shared_ptr<Grid> gridPtr_;
    Dune::GridPtr<Grid> dgfGridPtr_;

    std::shared_ptr<GridData> gridData_;
};

/*!
 * \brief The grid manager (this is the class used by the user) for all supported grid managers that constructs a grid
 *        from information in the input file and handles the data.
 * \note  This class is specialised below for all supported grid managers. It inherits the functionality of the base class.
 */
template <class Grid>
class GridManager
: public GridManagerBase<Grid>
{};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Specializations //////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
 * \brief Provides a grid manager for YaspGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.
 * The following keys are recognized:
 * - File : a DGF file to load the coarse grid from
 * - UpperRight : extension of the domain
 * - Cells : the number of cells in each direction
 * - Periodic : true or false for each direction
 * - Overlap : overlap size in cells
 * - Partitioning : a non-standard load-balancing, number of processors per direction
 * - KeepPyhsicalOverlap : whether to keep the physical overlap
 *     in physical size or in number of cells upon refinement
 * - Refinement : the number of global refines to apply initially.
 *
 */
template<class ct, int dim>
class GridManager<Dune::YaspGrid<dim, Dune::EquidistantCoordinates<ct, dim> >>
: public GridManagerBase<Dune::YaspGrid<dim, Dune::EquidistantCoordinates<ct, dim> > >
{
public:
    using Grid = typename Dune::YaspGrid<dim, Dune::EquidistantCoordinates<ct, dim> >;
    using ParentType = GridManagerBase<Grid>;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    void init(const std::string& modelParamGroup = "")
    {
        // First try to create it from a DGF file in GridParameterGroup.File
        if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromDgfFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"));
            postProcessing_(modelParamGroup);
            return;
        }

        // Then look for the necessary keys to construct from the input file
        else if (hasParamInGroup(modelParamGroup, "Grid.UpperRight"))
        {

            // get the upper right corner coordinates
            const auto upperRight = getParamFromGroup<Dune::FieldVector<ct, dim>>(modelParamGroup, "Grid.UpperRight");

            // number of cells in each direction
            std::array<int, dim> cells; cells.fill(1);
            cells = getParamFromGroup<std::array<int, dim>>(modelParamGroup, "Grid.Cells", cells);

            // \todo TODO periodic boundaries with yasp (the periodicity concept of yasp grid is currently not supported, use dune-spgrid)
            // const auto periodic = getParamFromGroup<std::bitset<dim>>(modelParamGroup, "Grid.Periodic", std::bitset<dim>());
            const std::bitset<dim> periodic;

            // get the overlap
            const int overlap =  getParamFromGroup<int>(modelParamGroup, "Grid.Overlap", 1);

            // make the grid
            if (!hasParamInGroup(modelParamGroup, "Grid.Partitioning"))
            {
                // construct using default load balancing
                ParentType::gridPtr() = std::make_shared<Grid>(upperRight, cells, periodic, overlap);
            }
            else
            {
                // construct using user defined partitioning
                const auto partitioning = getParamFromGroup<std::array<int, dim>>(modelParamGroup, "Grid.Partitioning");
                Dune::YaspFixedSizePartitioner<dim> lb(partitioning);
                ParentType::gridPtr() = std::make_shared<Grid>(upperRight, cells, periodic, overlap, typename Grid::CollectiveCommunicationType(), &lb);
            }
            postProcessing_(modelParamGroup);
        }

        // Didn't find a way to construct the grid
        else
        {
            const auto prefix = modelParamGroup == "" ? modelParamGroup : modelParamGroup + ".";
            DUNE_THROW(ParameterException, "Please supply one of the parameters "
                                           << prefix + "Grid.UpperRight"
                                           << ", or a grid file in " << prefix + "Grid.File");

        }
    }

private:
    /*!
     * \brief Postprocessing for YaspGrid
     */
    void postProcessing_(const std::string& modelParamGroup)
    {
        // Check if should refine the grid
        bool keepPhysicalOverlap = getParamFromGroup<bool>(modelParamGroup, "Grid.KeepPhysicalOverlap", true);
        ParentType::grid().refineOptions(keepPhysicalOverlap);
        ParentType::maybeRefineGrid(modelParamGroup);
        ParentType::loadBalance();
    }
};

/*!
 * \brief Provides a grid manager for YaspGrids with non-zero offset
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.
 * The following keys are recognized:
 * - LowerLeft : lower left corner coordinates
 * - UpperRight : upper right corner coordinates
 * - Cells : the number of cells in each direction
 * - Periodic : true or false for each direction
 * - Overlap : overlap size in cells
 * - Partitioning : a non-standard load-balancing, number of processors per direction
 * - KeepPyhsicalOverlap : whether to keep the physical overlap
 *     in physical size or in number of cells upon refinement
 * - Refinement : the number of global refines to apply initially.
 *
 */
template<class ct, int dim>
class GridManager<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<ct, dim>>>
: public GridManagerBase<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<ct, dim>>>
{
public:
    using Grid = typename Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<ct, dim>>;
    using ParentType = GridManagerBase<Grid>;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    void init(const std::string& modelParamGroup = "")
    {
        // First try to create it from a DGF file in GridParameterGroup.File
        if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromDgfFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"));
            postProcessing_(modelParamGroup);
            return;
        }

        // Then look for the necessary keys to construct from the input file
        else if (hasParamInGroup(modelParamGroup, "Grid.UpperRight"))
        {
            using GlobalPosition = Dune::FieldVector<ct, dim>;
            const auto upperRight = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.UpperRight");
            const auto lowerLeft = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.LowerLeft", GlobalPosition(0.0));

            // number of cells in each direction
            std::array<int, dim> cells; cells.fill(1);
            cells = getParamFromGroup<std::array<int, dim>>(modelParamGroup, "Grid.Cells", cells);

            // \todo TODO periodic boundaries with yasp (the periodicity concept of yasp grid is currently not supported, use dune-spgrid)
            // const auto periodic = getParamFromGroup<std::bitset<dim>>(modelParamGroup, "Grid.Periodic", std::bitset<dim>());
            const std::bitset<dim> periodic;

            // get the overlap dependent on some template parameters
            const int overlap = getParamFromGroup<int>(modelParamGroup, "Grid.Overlap", 1);

            // make the grid
            if (!hasParamInGroup(modelParamGroup, "Grid.Partitioning"))
            {
                // construct using default load balancing
                ParentType::gridPtr() = std::make_shared<Grid>(lowerLeft, upperRight, cells, periodic, overlap);
            }
            else
            {
                // construct using user defined partitioning
                const auto partitioning = getParamFromGroup<std::array<int, dim>>(modelParamGroup, "Grid.Partitioning");
                Dune::YaspFixedSizePartitioner<dim> lb(partitioning);
                ParentType::gridPtr() = std::make_shared<Grid>(lowerLeft, upperRight, cells, periodic, overlap, typename Grid::CollectiveCommunicationType(), &lb);
            }
            postProcessing_(modelParamGroup);
        }

        // Didn't find a way to construct the grid
        else
        {
            const auto prefix = modelParamGroup == "" ? modelParamGroup : modelParamGroup + ".";
            DUNE_THROW(ParameterException, "Please supply one of the parameters "
                                           << prefix + "Grid.UpperRight"
                                           << ", or a grid file in " << prefix + "Grid.File");

        }
    }

private:
    /*!
     * \brief Postprocessing for YaspGrid
     */
    void postProcessing_(const std::string& modelParamGroup)
    {
        // Check if should refine the grid
        const bool keepPhysicalOverlap = getParamFromGroup<bool>(modelParamGroup, "Grid.KeepPhysicalOverlap", true);
        ParentType::grid().refineOptions(keepPhysicalOverlap);
        ParentType::maybeRefineGrid(modelParamGroup);
        ParentType::loadBalance();
    }
};

/*!
 * \brief Provides a grid manager for YaspGrids with different zones and grading
 *
 * All keys are expected to be in group GridParameterGroup.
 * The following keys are recognized:
 * - Positions0 : position array for x-coordinate
 * - Positions1 : position array for y-coordinate
 * - Positions2 : position array for z-coordinate
 * - Cells0 : number of cells array for x-coordinate
 * - Cells1 : number of cells array for y-coordinate
 * - Cells2 : number of cells array for z-coordinate
 * - Grading0 : grading factor array for x-coordinate
 * - Grading1 : grading factor array for y-coordinate
 * - Grading2 : grading factor array for z-coordinate
 * - Verbosity : whether the grid construction should output to standard out
 * - Periodic : true or false for each direction
 * - Overlap : overlap size in cells
 * - Partitioning : a non-standard load-balancing, number of processors per direction
 * - KeepPyhsicalOverlap : whether to keep the physical overlap
 *     in physical size or in number of cells upon refinement
 * - Refinement : the number of global refines to apply initially.
 *
 * The grading factor \f$ g \f$ specifies the ratio between the next and the current cell size:
 * \f$ g = \frac{h_{i+1}}{h_i} \f$.
 * Negative grading factors are converted to
 * \f$ g = -\frac{1}{g_\textrm{negative}} \f$
 * to avoid issues with imprecise fraction numbers.
 */
template<class ctype, int dim>
class GridManager<Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ctype, dim> >>
: public GridManagerBase<Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ctype, dim> > >
{
public:
    using Grid = typename Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ctype, dim> >;
    using ParentType = GridManagerBase<Grid>;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    void init(const std::string& modelParamGroup = "")
    {
        // Only construction from the input file is possible
        // Look for the necessary keys to construct from the input file
        // The positions
        std::array<std::vector<ctype>, dim> positions;
        for (int i = 0; i < dim; ++i)
            positions[i] = getParamFromGroup<std::vector<ctype>>(modelParamGroup, "Grid.Positions" + std::to_string(i));

        // the number of cells (has a default)
        std::array<std::vector<int>, dim> cells;
        for (int i = 0; i < dim; ++i)
        {
            cells[i].resize(positions[i].size()-1, 1.0);
            cells[i] = getParamFromGroup<std::vector<int>>(modelParamGroup, "Grid.Cells" + std::to_string(i), cells[i]);
        }

        // grading factor (has a default)
        std::array<std::vector<ctype>, dim> grading;
        for (int i = 0; i < dim; ++i)
        {
            grading[i].resize(positions[i].size()-1, 1.0);
            grading[i] = getParamFromGroup<std::vector<ctype>>(modelParamGroup, "Grid.Grading" + std::to_string(i), grading[i]);
        }

        // call the generic function
        init(positions, cells, grading, modelParamGroup);
    }

    /*!
     * \brief Make the grid using input data not read from the input file.
     */
    void init(const std::array<std::vector<ctype>, dim>& positions,
              const std::array<std::vector<int>, dim>& cells,
              const std::array<std::vector<ctype>, dim>& grading,
              const std::string& modelParamGroup = "")
    {


        // Additional arameters (they have a default)
        const int overlap = getParamFromGroup<int>(modelParamGroup, "Grid.Overlap", 1);
        const bool verbose = getParamFromGroup<bool>(modelParamGroup, "Grid.Verbosity", false);
        // \todo TODO periodic boundaries with yasp (the periodicity concept of yasp grid is currently not supported, use dune-spgrid)
        // const auto periodic = getParamFromGroup<std::bitset<dim>>(modelParamGroup, "Grid.Periodic", std::bitset<dim>());
        const std::bitset<dim> periodic;

        // Some sanity checks
        for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            if (cells[dimIdx].size() + 1 != positions[dimIdx].size())
            {
                DUNE_THROW(Dune::RangeError, "Make sure to specify correct \"Cells\" and \"Positions\" arrays");
            }
            if (grading[dimIdx].size() + 1 != positions[dimIdx].size())
            {
                DUNE_THROW(Dune::RangeError, "Make sure to specify correct \"Grading\" and \"Positions\" arrays");
            }
            ctype temp = std::numeric_limits<ctype>::lowest();
            for (unsigned int posIdx = 0; posIdx < positions[dimIdx].size(); ++posIdx)
            {
                if (temp > positions[dimIdx][posIdx])
                {
                    DUNE_THROW(Dune::RangeError, "Make sure to specify a monotone increasing \"Positions\" array");
                }
                temp = positions[dimIdx][posIdx];
            }
        }

        const auto globalPositions = computeGlobalPositions_(positions, cells, grading, verbose);

        // make the grid
        if (!hasParamInGroup(modelParamGroup, "Grid.Partitioning"))
        {
            // construct using default load balancing
            ParentType::gridPtr() = std::make_shared<Grid>(globalPositions, periodic, overlap);
        }
        else
        {
            // construct using user defined partitioning
            const auto partitioning = getParamFromGroup<std::array<int, dim>>(modelParamGroup, "Grid.Partitioning");
            Dune::YaspFixedSizePartitioner<dim> lb(partitioning);
            ParentType::gridPtr() = std::make_shared<Grid>(globalPositions, periodic, overlap, typename Grid::CollectiveCommunicationType(), &lb);
        }

        postProcessing_(modelParamGroup);
    }

private:
    /*!
     * \brief Postprocessing for YaspGrid
     */
    void postProcessing_(const std::string& modelParamGroup)
    {
        // Check if should refine the grid
        const bool keepPhysicalOverlap = getParamFromGroup<bool>(modelParamGroup, "Grid.KeepPhysicalOverlap", true);
        ParentType::grid().refineOptions(keepPhysicalOverlap);
        ParentType::maybeRefineGrid(modelParamGroup);
        ParentType::loadBalance();
    }

    //! Compute the global position tensor grid from the given positions, cells, and grading factors
    std::array<std::vector<ctype>, dim>
    computeGlobalPositions_(const std::array<std::vector<ctype>, dim>& positions,
                            const std::array<std::vector<int>, dim>& cells,
                            const std::array<std::vector<ctype>, dim>& grading,
                            bool verbose = false)
    {
        std::array<std::vector<ctype>, dim> globalPositions;
        using std::pow;
        for (int dimIdx = 0; dimIdx < dim; dimIdx++)
        {
            for (int zoneIdx = 0; zoneIdx < cells[dimIdx].size(); ++zoneIdx)
            {
                ctype lower = positions[dimIdx][zoneIdx];
                ctype upper = positions[dimIdx][zoneIdx+1];
                int numCells = cells[dimIdx][zoneIdx];
                ctype gradingFactor = grading[dimIdx][zoneIdx];
                ctype length = upper - lower;
                ctype height = 1.0;
                bool increasingCellSize = false;

                if (verbose)
                {
                    std::cout << "dim " << dimIdx
                              << " lower "  << lower
                              << " upper "  << upper
                              << " numCells "  << numCells
                              << " grading "  << gradingFactor;
                }

                if (gradingFactor > 1.0)
                {
                    increasingCellSize = true;
                }

                // take absolute values and reverse cell size increment to achieve
                // reverse behavior for negative values
                if (gradingFactor < 0.0)
                {
                    using std::abs;
                    gradingFactor = abs(gradingFactor);
                    if (gradingFactor < 1.0)
                    {
                        increasingCellSize = true;
                    }
                }

                // if the grading factor is exactly 1.0 do equal spacing
                if (gradingFactor > 1.0 - 1e-7 && gradingFactor < 1.0 + 1e-7)
                {
                    height = 1.0 / numCells;
                    if (verbose)
                    {
                        std::cout << " -> h "  << height * length << std::endl;
                    }
                }
                // if grading factor is not 1.0, do power law spacing
                else
                {
                    height = (1.0 - gradingFactor) / (1.0 - pow(gradingFactor, numCells));

                    if (verbose)
                    {
                        std::cout << " -> grading_eff "  << gradingFactor
                                  << " h_min "  << height * pow(gradingFactor, 0) * length
                                  << " h_max "  << height * pow(gradingFactor, numCells-1) * length
                                  << std::endl;
                    }
                }

                std::vector<ctype> localPositions;
                localPositions.push_back(0);
                for (int i = 0; i < numCells-1; i++)
                {
                    ctype hI = height;
                    if (!(gradingFactor < 1.0 + 1e-7 && gradingFactor > 1.0 - 1e-7))
                    {
                        if (increasingCellSize)
                        {
                            hI *= pow(gradingFactor, i);
                        }
                        else
                        {
                            hI *= pow(gradingFactor, numCells-i-1);
                        }
                    }
                    localPositions.push_back(localPositions[i] + hI);
                }

                for (int i = 0; i < localPositions.size(); i++)
                {
                    localPositions[i] *= length;
                    localPositions[i] += lower;
                }

                for (unsigned int i = 0; i < localPositions.size(); ++i)
                {
                    globalPositions[dimIdx].push_back(localPositions[i]);
                }
            }
            globalPositions[dimIdx].push_back(positions[dimIdx].back());
        }

        return globalPositions;
    }
};

/*!
 * \brief Provides a grid manager for OneDGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.
 * The following keys are recognized:
 * - LeftBoundary : start coordinate
 * - RightBoundary : end coordinate
 * - Cells : the number of cell
 * - RefinementType : local or copy
 * - Refinement : the number of global refines to apply initially.
 *
 */
template<>
class GridManager<Dune::OneDGrid>
: public GridManagerBase<Dune::OneDGrid>
{
public:
    using Grid = Dune::OneDGrid;
    using ParentType = GridManagerBase<Grid>;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    void init(const std::string& modelParamGroup = "")
    {

        // try to create it from a DGF or msh file in GridParameterGroup.File
        if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromDgfFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"));
            postProcessing_(modelParamGroup);
            return;
        }

        // Look for the necessary keys to construct from the input file
        else if (hasParamInGroup(modelParamGroup, "Grid.RightBoundary"))
        {
            // The required parameters
            using CoordinateType = typename Grid::ctype;
            const auto leftBoundary = getParamFromGroup<CoordinateType>(modelParamGroup, "Grid.LeftBoundary", 0.0);
            const auto rightBoundary = getParamFromGroup<CoordinateType>(modelParamGroup, "Grid.RightBoundary");
            const int cells = getParamFromGroup<int>(modelParamGroup, "Grid.Cells", 1);

            ParentType::gridPtr() = std::make_shared<Grid>(cells, leftBoundary, rightBoundary);
            postProcessing_(modelParamGroup);
            return;
        }

        // Look for the necessary keys to construct from the input file with just a coordinates vector
        else if (hasParamInGroup(modelParamGroup, "Grid.Coordinates"))
        {
            const auto coordinates = getParamFromGroup<std::vector<typename Grid::ctype>>(modelParamGroup, "Grid.Coordinates");
            ParentType::gridPtr() = std::make_shared<Grid>(coordinates);
            postProcessing_(modelParamGroup);
        }

        // Didn't find a way to construct the grid
        else
        {
            const auto prefix = modelParamGroup == "" ? modelParamGroup : modelParamGroup + ".";
            DUNE_THROW(ParameterException, "Please supply one of the parameters "
                                           << prefix + "Grid.RightBoundary"
                                           << ", or " << prefix + "Grid.Coordinates"
                                           << ", or a grid file in " << prefix + "Grid.File");
        }
    }

    /*!
     * \brief Call loadBalance() function of the grid. OneDGrid is not parallel an thus cannot communicate.
     */
    void loadBalance() {}

private:
    /*!
     * \brief Do some operatrion after making the grid, like global refinement
     */
    void postProcessing_(const std::string& modelParamGroup)
    {
        // Set refinement type
        const auto refType = getParamFromGroup<std::string>(modelParamGroup, "Grid.RefinementType", "Local");
        if (refType == "Local")
            ParentType::grid().setRefinementType(Dune::OneDGrid::RefinementType::LOCAL);
        else if (refType == "Copy")
            ParentType::grid().setRefinementType(Dune::OneDGrid::RefinementType::COPY);
        else
            DUNE_THROW(Dune::IOError, "OneGrid only supports 'Local' or 'Copy' as refinment type. Not '"<< refType<<"'!");

        // Check if should refine the grid
        ParentType::maybeRefineGrid(modelParamGroup);
        loadBalance();
    }
};

#if HAVE_UG

/*!
 * \brief Provides a grid manager for UGGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.

 * The following keys are recognized:
 * - File : A DGF or gmsh file to load from, type detection by file extension
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 * - CellType : "Cube" or "Simplex" to be used for structured grids
 * - Refinement : the number of global refines to perform
 * - Verbosity : whether the grid construction should output to standard out
 * - HeapSize: The heapsize used to allocate memory
 * - BoundarySegments : whether to insert boundary segments into the grid
 *
 */
template<int dim>
class GridManager<Dune::UGGrid<dim>>
: public GridManagerBase<Dune::UGGrid<dim>>
{
public:
    using Grid = typename Dune::UGGrid<dim>;
    using ParentType = GridManagerBase<Grid>;
    using Element = typename Grid::template Codim<0>::Entity;

    /*!
     * \brief Make the UGGrid.
     */
    void init(const std::string& modelParamGroup = "")
    {

        // try to create it from a DGF or msh file in GridParameterGroup.File
        if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            preProcessing_(modelParamGroup);
            ParentType::makeGridFromFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"), modelParamGroup);
            postProcessing_(modelParamGroup);
            return;
        }

        // Then look for the necessary keys to construct from the input file
        else if (hasParamInGroup(modelParamGroup, "Grid.UpperRight"))
        {
            preProcessing_(modelParamGroup);
            // make the grid
            const auto cellType = getParamFromGroup<std::string>(modelParamGroup, "Grid.CellType", "Cube");
            if (cellType == "Cube")
                ParentType::template makeStructuredGrid<dim, dim>(ParentType::CellType::Cube, modelParamGroup);
            else if (cellType == "Simplex")
                ParentType::template makeStructuredGrid<dim, dim>(ParentType::CellType::Simplex, modelParamGroup);
            else
                DUNE_THROW(Dune::IOError, "UGGrid only supports 'Cube' or 'Simplex' as cell type. Not '"<< cellType<<"'!");
            postProcessing_(modelParamGroup);
        }

        // Didn't find a way to construct the grid
        else
        {
            const auto prefix = modelParamGroup == "" ? modelParamGroup : modelParamGroup + ".";
            DUNE_THROW(ParameterException, "Please supply one of the parameters "
                                           << prefix + "Grid.UpperRight"
                                           << ", or a grid file in " << prefix + "Grid.File");

        }
    }

    /*!
     * \brief Call loadBalance() function of the grid.
     * \note This overload is necessary because UGGrid cannot load balance element data yet.
     *       For dgf grids using element data in parallel will (hopefully) through an error
     *       For gmsh grids the parameters are read on every process so we use too much memory but
     *       the parameters are available via the insertionIndex of the level 0 element.
     */
    void loadBalance()
    {
        if (Dune::MPIHelper::getCollectiveCommunication().size() > 1)
        {
            // if we may have dgf parameters use load balancing of the dgf pointer
            if(ParentType::enableDgfGridPointer_)
            {
                ParentType::dgfGridPtr().loadBalance();
                // update the grid data
                ParentType::gridData_ = std::make_shared<typename ParentType::GridData>(ParentType::dgfGridPtr());
            }

            // if we have gmsh parameters we have to manually load balance the data
            else if (ParentType::enableGmshDomainMarkers_)
            {
                // element and face markers are communicated during load balance
                auto dh = ParentType::gridData_->createGmshDataHandle();
                ParentType::gridPtr()->loadBalance(dh.interface());
                // Right now, UGGrid cannot communicate element data. If this gets implemented, communicate the data here:
                // ParentType::gridPtr()->communicate(dh.interface(), Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
            }
            else
                ParentType::gridPtr()->loadBalance();
        }
    }

private:
    /*!
     * \brief Do some operations before making the grid
     */
    void preProcessing_(const std::string& modelParamGroup)
    {
        if(hasParamInGroup(modelParamGroup, "Grid.HeapSize"))
            Grid::setDefaultHeapSize(getParamFromGroup<unsigned>(modelParamGroup, "Grid.HeapSize"));
    }

    /*!
     * \brief Do some operatrion after making the grid, like global refinement
     */
    void postProcessing_(const std::string& modelParamGroup)
    {
        // Set refinement type
        const auto refType = getParamFromGroup<std::string>(modelParamGroup, "Grid.RefinementType", "Local");
        if (refType == "Local")
            ParentType::grid().setRefinementType(Dune::UGGrid<dim>::RefinementType::LOCAL);
        else if (refType == "Copy")
            ParentType::grid().setRefinementType(Dune::UGGrid<dim>::RefinementType::COPY);
        else
            DUNE_THROW(Dune::IOError, "UGGrid only supports 'Local' or 'Copy' as refinment type. Not '"<< refType<<"'!");

        // Set closure type
        const auto closureType = getParamFromGroup<std::string>(modelParamGroup, "Grid.ClosureType", "Green");
        if (closureType == "Green")
            ParentType::grid().setClosureType(Dune::UGGrid<dim>::ClosureType::GREEN);
        else if (closureType == "None")
            ParentType::grid().setClosureType(Dune::UGGrid<dim>::ClosureType::NONE);
        else
            DUNE_THROW(Dune::IOError, "UGGrid only supports 'Green' or 'None' as closure type. Not '"<< closureType<<"'!");

        // Check if should refine the grid
        ParentType::maybeRefineGrid(modelParamGroup);
        // do load balancing
        loadBalance();
    }
};

#endif // HAVE_UG
//
#if HAVE_DUNE_ALUGRID

/*!
 * \brief Provides a grid manager for Dune ALUGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.

 * The following keys are recognized:
 * - File : A DGF or gmsh file to load from, type detection by file extension
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 * - Refinement : the number of global refines to perform
 * - Verbosity : whether the grid construction should output to standard out
 * - BoundarySegments : whether to insert boundary segments into the grid
 *
 */
template<int dim, int dimworld, Dune::ALUGridElementType elType, Dune::ALUGridRefinementType refinementType>
class GridManager<Dune::ALUGrid<dim, dimworld, elType, refinementType>>
: public GridManagerBase<Dune::ALUGrid<dim, dimworld, elType, refinementType>>
{
public:
    using Grid = Dune::ALUGrid<dim, dimworld, elType, refinementType>;
    using ParentType = GridManagerBase<Grid>;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    void init(const std::string& modelParamGroup = "", bool adaptiveRestart = false)
    {
        // restarting an adaptive grid using Dune's BackupRestoreFacility
        // TODO: the part after first || is backward compatibilty with old sequential models remove once sequential adpative restart is replaced
        if (adaptiveRestart || hasParam("Restart") || hasParam("TimeManager.Restart"))
        {
            auto restartTime = getParamFromGroup<double>(modelParamGroup, "TimeLoop.Restart", 0.0);
            // TODO: backward compatibilty with old sequential models remove once sequential adpative restart is replaced
            if (hasParam("Restart") || hasParam("TimeManager.Restart"))
            {
                restartTime = getParamFromGroup<double>("TimeManager", "Restart");
                std::cerr << "Warning: You are using a deprecated restart mechanism. The usage will change in the future.\n";
            }

            const int rank = Dune::MPIHelper::getCollectiveCommunication().rank();
            const std::string name = getParamFromGroup<std::string>(modelParamGroup, "Problem.Name");
            std::ostringstream oss;
            oss << name << "_time=" << restartTime << "_rank=" << rank << ".grs";
            std::cout << "Restoring an ALUGrid from " << oss.str() << std::endl;
            ParentType::gridPtr() = std::shared_ptr<Grid>(Dune::BackupRestoreFacility<Grid>::restore(oss.str()));
            ParentType::loadBalance();
            return;
        }

        // try to create it from a DGF or msh file in GridParameterGroup.File
        else if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            makeGridFromFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"), modelParamGroup);
            ParentType::maybeRefineGrid(modelParamGroup);
            ParentType::loadBalance();
            return;
        }

        // Then look for the necessary keys to construct from the input file
        else if (hasParamInGroup(modelParamGroup, "Grid.UpperRight"))
        {
            // make a structured grid
            if (elType == Dune::cube)
                ParentType::template makeStructuredGrid<dim, dimworld>(ParentType::CellType::Cube, modelParamGroup);
            else if (elType == Dune::simplex)
                ParentType::template makeStructuredGrid<dim, dimworld>(ParentType::CellType::Simplex, modelParamGroup);
            else
                DUNE_THROW(Dune::IOError, "ALUGrid only supports Dune::cube or Dune::simplex as cell type!");
            ParentType::maybeRefineGrid(modelParamGroup);
            ParentType::loadBalance();
        }

        // Didn't find a way to construct the grid
        else
        {
            const auto prefix = modelParamGroup == "" ? modelParamGroup : modelParamGroup + ".";
            DUNE_THROW(ParameterException, "Please supply one of the parameters "
                                           << prefix + "Grid.UpperRight"
                                           << ", or a grid file in " << prefix + "Grid.File");

        }
    }

    /*!
     * \brief Makes a grid from a file. We currently support *.dgf (Dune Grid Format) and *.msh (Gmsh mesh format).
     */
    void makeGridFromFile(const std::string& fileName,
                          const std::string& modelParamGroup)
    {
        // We found a file in the input file...does it have a supported extension?
        const std::string extension = ParentType::getFileExtension(fileName);
        if(extension != "dgf" && extension != "msh")
            DUNE_THROW(Dune::IOError, "Grid type " << Dune::className<Grid>() << " only supports DGF (*.dgf) and Gmsh (*.msh) grid files but the specified filename has extension: *."<< extension);

        // make the grid
        if (extension == "dgf")
        {
            ParentType::enableDgfGridPointer_ = true;
            ParentType::dgfGridPtr() = Dune::GridPtr<Grid>(fileName.c_str(), Dune::MPIHelper::getCommunicator());
            ParentType::gridData_ = std::make_shared<typename ParentType::GridData>(ParentType::dgfGridPtr());
        }
        if (extension == "msh")
        {
            // get some optional parameters
            const bool verbose = getParamFromGroup<bool>(modelParamGroup, "Grid.Verbosity", false);
            const bool boundarySegments = getParamFromGroup<bool>(modelParamGroup, "Grid.BoundarySegments", false);
            const bool domainMarkers = getParamFromGroup<bool>(modelParamGroup, "Grid.DomainMarkers", false);

            if (domainMarkers)
                ParentType::enableGmshDomainMarkers_ = true;

            // only filll the factory for rank 0
            if (domainMarkers)
            {
                std::vector<int> boundaryMarkersInsertionIndex, boundaryMarkers, faceMarkers, elementMarkers;
                auto gridFactory = std::make_unique<Dune::GridFactory<Grid>>();
                if (Dune::MPIHelper::getCollectiveCommunication().rank() == 0)
                    Dune::GmshReader<Grid>::read(*gridFactory, fileName, boundaryMarkersInsertionIndex, elementMarkers, verbose, boundarySegments);

                ParentType::gridPtr() = std::shared_ptr<Grid>(gridFactory->createGrid());

                // reorder boundary markers according to boundarySegmentIndex
                boundaryMarkers.resize(ParentType::gridPtr()->numBoundarySegments(), 0);
                faceMarkers.resize(ParentType::gridPtr()->leafGridView().size(1), 0);
                const auto& indexSet = ParentType::gridPtr()->leafGridView().indexSet();
                for (const auto& element : elements(ParentType::gridPtr()->leafGridView()))
                {
                    for (const auto& intersection : intersections(ParentType::gridPtr()->leafGridView(), element))
                    {
                        if (intersection.boundary() && gridFactory->wasInserted(intersection))
                        {
                            auto marker = boundaryMarkersInsertionIndex[gridFactory->insertionIndex(intersection)];
                            boundaryMarkers[intersection.boundarySegmentIndex()] = marker;
                            faceMarkers[indexSet.index(element.template subEntity<1>(intersection.indexInInside()))] = marker;
                        }
                    }
                }

                ParentType::gridData_ = std::make_shared<typename ParentType::GridData>(ParentType::gridPtr(), std::move(gridFactory),
                                                       std::move(elementMarkers), std::move(boundaryMarkers), std::move(faceMarkers));
            }
            else
            {
                auto gridFactory = std::make_unique<Dune::GridFactory<Grid>>();
                if (Dune::MPIHelper::getCollectiveCommunication().rank() == 0)
                    Dune::GmshReader<Grid>::read(*gridFactory, fileName, verbose, boundarySegments);

                ParentType::gridPtr() = std::shared_ptr<Grid>(gridFactory->createGrid());
            }
        }
    }
};

#endif // HAVE_DUNE_ALUGRID

#if HAVE_DUNE_FOAMGRID

/*!
 * \brief Provides a grid manager for FoamGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.

 * The following keys are recognized:
 * - File : A DGF or gmsh file to load from, type detection by file extension
 * - Verbosity : whether the grid construction should output to standard out
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 *
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
        // try to create it from file
        if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"), modelParamGroup);
            ParentType::maybeRefineGrid(modelParamGroup);
            ParentType::loadBalance();
            return;
        }

        // Then look for the necessary keys to construct a structured grid from the input file
        else if (hasParamInGroup(modelParamGroup, "Grid.UpperRight"))
        {
            ParentType::template makeStructuredGrid<dim, dimworld>(ParentType::CellType::Simplex, modelParamGroup);
            ParentType::maybeRefineGrid(modelParamGroup);
            ParentType::loadBalance();
        }

        // Didn't find a way to construct the grid
        else
        {
            const auto prefix = modelParamGroup == "" ? modelParamGroup : modelParamGroup + ".";
            DUNE_THROW(ParameterException, "Please supply one of the parameters "
                                           << prefix + "Grid.UpperRight"
                                           << ", or a grid file in " << prefix + "Grid.File");

        }
    }
};

/*!
 * \brief Provides a grid manager for FoamGrids of dim 1
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.

 * The following keys are recognized:
 * - File : A DGF or gmsh file to load from, type detection by file extension
 * - Verbosity : whether the grid construction should output to standard out
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 *
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
        // try to create it from file
        if (hasParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"), modelParamGroup);
            ParentType::maybeRefineGrid(modelParamGroup);
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
        ParentType::loadBalance();
    }
};

#endif // HAVE_DUNE_FOAMGRID

#if HAVE_DUNE_SPGRID

/*!
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
            const auto overlap = getParamFromGroup<int>(paramGroup, "Grid.Overlap", 1);
            IntArray spOverlap; spOverlap.fill(overlap);

            using Domain = typename Grid::Domain;
            std::vector< typename Domain::Cube > cubes;
            cubes.push_back( typename Domain::Cube( lowerLeft, upperRight ) );
            Domain domain( cubes, typename Domain::Topology( static_cast<unsigned int>(periodic.to_ulong()) ) );
            ParentType::gridPtr() = std::make_shared<Grid>( domain, cells, spOverlap );
            ParentType::maybeRefineGrid(paramGroup);
            ParentType::loadBalance();
        }
        else
        {
            const auto prefix = paramGroup == "" ? paramGroup : paramGroup + ".";
            DUNE_THROW(ParameterException, "Please supply a grid file in " << prefix << "Grid.File or " << prefix << "Grid.UpperRight/Cells.");
        }
    }
};

#endif // HAVE_DUNE_SPGRID

} // end namespace Dumux

#endif
