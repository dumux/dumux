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
 * \brief Provides a grid manager for pore-network models
 */
#ifndef DUMUX_PORE_NETWORK_GRID_MANAGER_HH
#define DUMUX_PORE_NETWORK_GRID_MANAGER_HH

#include <iostream>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>

// FoamGrid specific includes
#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/dgffoam.hh>
#else
static_assert(false, "dune-foamgrid required!");
#endif

#include <dumux/common/parameters.hh>

#include "griddata.hh"
#include "structuredlatticegridcreator.hh"

#include "dgfwriter.hh"

namespace Dumux {

/*!
 * \ingroup PoreNetworkModel
 * \brief A grid manager for pore-network models
 */
template<int dimWorld, class GData = Dumux::PoreNetworkGridData<Dune::FoamGrid<1, dimWorld>>>
class PoreNetworkGridManager
{
    using GridType = Dune::FoamGrid<1, dimWorld>;
    using Element = typename GridType::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using Grid = GridType;
    using GridData = GData;

    void init(const std::string& paramGroup = "")
    {
        Dune::Timer timer;
        paramGroup_ = paramGroup;

        // First try to create it from a DGF file in GridParameterGroup.File
        if (hasParamInGroup(paramGroup, "Grid.File"))
        {
            makeGridFromDgfFile(getParamFromGroup<std::string>(paramGroup, "Grid.File"));
            loadBalance();
        }
        else // no grid file found
        {
            try
            {
                createOneDGrid_();
            }
            catch(Dune::RangeError& e)
            {
                makeGridFromStructuredLattice();
            }

            loadBalance();
        }
        // write the final grid to a dgf file, if desired
        if (getParamFromGroup<bool>(paramGroup_, "Grid.WriteDgfFile", false))
        {
            const auto defaultName = enableDgfGridPointer_ ? "pnm-grid-from-dgf.dgf" : "pnm-grid-from-factory.dgf";
            const auto fileName = getParamFromGroup<std::string>(paramGroup, "Grid.DgfName", defaultName);
            writeDgf(fileName, grid().leafGridView(), *gridData_);
        }

        timer.stop();
        const auto gridType = enableDgfGridPointer_ ?  "grid from dgf file" : "generic grid from structured lattice";
        std::cout << "Creating "  << gridType << " with " << grid().leafGridView().size(0) << " elements and "
                  << grid().leafGridView().size(1) << " vertices took " << timer.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief Make a grid from a DGF file.
     */
    void makeGridFromDgfFile(const std::string& fileName)
    {
        // We found a file in the input file...does it have a supported extension?
        const std::string extension = getFileExtension(fileName);
        if (extension != "dgf")
            DUNE_THROW(Dune::IOError, "Grid type " << Dune::className<Grid>() << " only supports DGF (*.dgf) but the specified filename has extension: *."<< extension);

        enableDgfGridPointer_ = true;
        dgfGridPtr() = Dune::GridPtr<Grid>(fileName.c_str(), Dune::MPIHelper::getCommunicator());
        gridData_ = std::make_shared<GridData>(dgfGridPtr_, paramGroup_);

        if (getParamFromGroup<bool>(paramGroup_, "Grid.Sanitize", false))
            sanitizeGrid();
    }

    /*!
     * \brief Make a grid based on a structured lattice which allows to randomly delete elements based on Raoof et al. 2009
     */
    void makeGridFromStructuredLattice()
    {
        StructuredLatticeGridCreator<dimWorld> creator;
        creator.init(paramGroup_);
        gridPtr() = creator.gridPtr();

        gridData_ = std::make_shared<GridData>(gridPtr_, paramGroup_);

        if (getParamFromGroup<bool>(paramGroup_, "Grid.Sanitize", true))
            sanitizeGrid();
        else
            std::cout << "\nWARNING: Set Grid.Sanitize = true in order to remove insular patches of elements not connected to the boundary." << std::endl;

        gridData_->assignParameters();
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
     * \brief Returns a reference to the grid.
     */
    Grid& grid()
    { return enableDgfGridPointer_ ? *dgfGridPtr() : *gridPtr(); }

    /*!
     * \brief Returns the grid data.
     */
    std::shared_ptr<GridData> getGridData() const
    {
        if (!gridData_)
            DUNE_THROW(Dune::IOError, "No grid data available");

        return gridData_;
    }

    /*!
     * \brief Call loadBalance() function of the grid.
     */
    void loadBalance()
    {
        if (Dune::MPIHelper::getCollectiveCommunication().size() > 1)
        {
            // if we may have dgf parameters use load balancing of the dgf pointer
            if (enableDgfGridPointer_)
            {
                dgfGridPtr().loadBalance();
                // update the grid data
                gridData_ = std::make_shared<GridData>(dgfGridPtr(), paramGroup_);
            }

            else
                gridPtr()->loadBalance();
        }
    }

    /*!
     * \brief post-processing to remove insular groups of elements that are not connected to a Dirichlet boundary
     */
    void sanitizeGrid()
    {
        // evaluate the coordination numbers to check if all pores are connected to at least one throat
        gridData_->getCoordinationNumbers();

        // for dgf grid, copy the data to peristent containers
        if (enableDgfGridPointer_)
            gridData_->copyDgfData();

        // pruning -- remove elements not connected to a Dirichlet boundary (marker == 1)
        const auto pruningSeedIndices = getParamFromGroup<std::vector<int>>(paramGroup_, "Grid.PruningSeedIndices", std::vector<int>{1});
        const auto gridView = grid().leafGridView();
        std::vector<bool> elementIsConnected(gridView.size(0), false);

        auto boundaryIdx = [&](const auto& vertex)
        {
            if (enableDgfGridPointer_)
                return static_cast<int>(dgfGridPtr_.parameters(vertex)[gridData_->parameterIndex("PoreLabel")]);
            else
                return static_cast<int>(gridData_->poreLabelAtPosForGenericGrid(vertex.geometry().center()));
        };

        for (const auto& element : elements(gridView))
        {
            const auto eIdx = gridView.indexSet().index(element);
            if (elementIsConnected[eIdx])
                continue;

            // try to find a seed from which to start the search process
            bool isSeed = false;
            bool hasNeighbor = false;
            for (const auto& intersection : intersections(gridView, element))
            {
                auto vertex = element.template subEntity<1>(intersection.indexInInside());
                // seed found
                if (std::any_of(pruningSeedIndices.begin(), pruningSeedIndices.end(),
                               [&]( const int i ){ return i == boundaryIdx(vertex); }))
                {
                    isSeed = true;
                    // break;
                }

                if (intersection.neighbor())
                    hasNeighbor = true;
            }

            if (!hasNeighbor)
                continue;

            if (isSeed)
            {
                elementIsConnected[eIdx] = true;

                // use iteration instead of recursion here because the recursion can get too deep
                std::stack<Element> elementStack;
                elementStack.push(element);
                while (!elementStack.empty())
                {
                    auto e = elementStack.top();
                    elementStack.pop();
                    for (const auto& intersection : intersections(gridView, e))
                    {
                        if (!intersection.boundary())
                        {
                            auto outside = intersection.outside();
                            auto nIdx = gridView.indexSet().index(outside);
                            if (!elementIsConnected[nIdx])
                            {
                                elementIsConnected[nIdx] = true;
                                elementStack.push(outside);
                            }
                        }
                    }
                }
            }
        }

        if (std::none_of(elementIsConnected.begin(), elementIsConnected.end(), [](const bool i){return i;}))
            DUNE_THROW(Dune::InvalidStateException, "sanitize() deleted all elements! Check your boundary face indices. "
                << "If the problem persisits, add at least one of your boundary face indices to PruningSeedIndices");

        // remove the elements in the grid
        std::size_t numDeletedElements = 0;
        grid().preGrow();
        for (const auto& element : elements(gridView))
        {
            const auto eIdx = gridView.indexSet().index(element);
            if (!elementIsConnected[eIdx])
            {
                grid().removeElement(element);
                ++numDeletedElements;
            }
        }
        // triggers the grid growth process
        grid().grow();
        grid().postGrow();

        // resize the parameters for dgf grids
        if (enableDgfGridPointer_)
            gridData_->resizeParameterContainers();

        if (numDeletedElements > 0)
            std::cout << "\nDeleted " << numDeletedElements << " elements not connected to a boundary." << std::endl;
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
    * \brief A state variable if the DGF Dune::GridPtr has been enabled.
    *        It is always enabled if a DGF grid file was used to create the grid.
    */
    bool enableDgfGridPointer_ = false;

    std::shared_ptr<Grid> gridPtr_;
    Dune::GridPtr<Grid> dgfGridPtr_;

    std::shared_ptr<GridData> gridData_;

    std::string paramGroup_;

private:
    void createOneDGrid_()
    {
        const auto lowerLeft = getParamFromGroup<GlobalPosition>(paramGroup_, "Grid.LowerLeft", GlobalPosition(0.0));
        const auto upperRight = getParamFromGroup<GlobalPosition>(paramGroup_, "Grid.UpperRight");
        const auto numPores = getParamFromGroup<unsigned int>(paramGroup_, "Grid.NumPores");
        const auto cells = numPores - 1;

        // create a step vector
        GlobalPosition step = upperRight;
        step -= lowerLeft, step /= cells;

        // make the grid (structured interval grid in dimworld space)
        Dune::GridFactory<Grid> factory;

        // create the vertices
        GlobalPosition globalPos = lowerLeft;
        for (unsigned int vIdx = 0; vIdx <= cells; vIdx++, globalPos += step)
            factory.insertVertex(globalPos);

        // create the cells
        for(unsigned int eIdx = 0; eIdx < cells; eIdx++)
            factory.insertElement(Dune::GeometryTypes::line, {eIdx, eIdx+1});

        gridPtr() = std::shared_ptr<Grid>(factory.createGrid());
        gridData_ = std::make_shared<GridData>(gridPtr_, paramGroup_);
        gridData_->assignParameters();
    }
};

}

#endif
