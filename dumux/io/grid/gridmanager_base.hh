// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 */
#ifndef DUMUX_IO_GRID_MANAGER_BASE_HH
#define DUMUX_IO_GRID_MANAGER_BASE_HH

#include <array>
#include <bitset>
#include <memory>
#include <sstream>

#include <dune/common/exceptions.hh>
#include <dune/common/classname.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/io/vtk/vtkreader.hh>

#include "griddata.hh"

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief The grid manager (this is the class used by the user) for all supported grid managers that constructs a grid
 *        from information in the input file and handles the data.
 * \note  This class is specialised below for all supported grid managers. It inherits the functionality of the base class Dumux::GridManagerBase.
 */
template <class Grid>
class GridManager;

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
                    "The header with the GridManager specialization for grid type " << Dune::className<Grid>()
                    << " is not included or no specialization has been implemented!"
                    << " In case of the latter, consider providing your own GridManager.");
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
        if (!gridData_)
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
     * \brief Makes a grid from a file. We currently support
     *     - dgf (Dune Grid Format)
     *     - msh (Gmsh mesh format)
     *     - vtp/vtu (VTK file formats)
     */
    void makeGridFromFile(const std::string& fileName,
                          const std::string& modelParamGroup)
    {
        // We found a file in the input file...does it have a supported extension?
        const std::string extension = getFileExtension(fileName);
        if (extension != "dgf" && extension != "msh" && extension != "vtu" && extension != "vtp")
            DUNE_THROW(Dune::IOError, "Grid type " << Dune::className<Grid>() << " doesn't support grid files with extension: *."<< extension);

        // Dune Grid Format (DGF) files
        if (extension == "dgf")
        {
            enableDgfGridPointer_ = true;
            dgfGridPtr() = Dune::GridPtr<Grid>(fileName.c_str(), Dune::MPIHelper::getCommunicator());
            gridData_ = std::make_shared<GridData>(dgfGridPtr_);
        }

        // Gmsh mesh format
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

        // VTK file formats for unstructured grids
        else if (extension == "vtu" || extension == "vtp")
        {
            if (Dune::MPIHelper::getCollectiveCommunication().size() > 1)
                DUNE_THROW(Dune::NotImplemented, "Reading grids in parallel from VTK file formats is currently not supported!");

            VTKReader vtkReader(fileName);
            VTKReader::Data cellData, pointData;
            auto gridFactory = std::make_unique<Dune::GridFactory<Grid>>();
            const bool verbose = getParamFromGroup<bool>(modelParamGroup, "Grid.Verbosity", false);
            gridPtr() = vtkReader.readGrid(*gridFactory, cellData, pointData, verbose);
            gridData_ = std::make_shared<GridData>(gridPtr_, std::move(gridFactory), std::move(cellData), std::move(pointData));
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

template <class Grid>
class GridManager : public GridManagerBase<Grid> {};

} // end namespace Dumux

#endif
