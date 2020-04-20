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
 * \brief Grid manager specialization for ALUGrid
 */
#ifndef DUMUX_IO_GRID_MANAGER_ALU_HH
#define DUMUX_IO_GRID_MANAGER_ALU_HH

// ALUGrid specific includes
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#endif

#ifndef DUMUX_IO_GRID_MANAGER_BASE_HH
#include <dumux/io/grid/gridmanager_base.hh>
#endif

#include <dumux/common/boundaryflag.hh>

namespace Dumux {

#if HAVE_DUNE_ALUGRID

/*!
 * \ingroup InputOutput
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
            if (elType == Dune::cube)
                makeStructuredGrid<dim, dimworld>(ParentType::CellType::Cube, modelParamGroup);
            else if (elType == Dune::simplex)
                makeStructuredGrid<dim, dimworld>(ParentType::CellType::Simplex, modelParamGroup);
            else
                DUNE_THROW(Dune::IOError, "ALUGrid only supports Dune::cube or Dune::simplex as cell type!");

            ParentType::maybeRefineGrid(modelParamGroup);
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

    /*!
     * \brief Makes a structured cube grid using the structured grid factory
     */
    template <int dimension, int dimensionworld, std::enable_if_t<dimension != dimensionworld, int> = 0>
    void makeStructuredGrid(typename ParentType::CellType cellType,
                            const std::string& modelParamGroup)
    {
        DUNE_THROW(Dune::IOError, "ALUGrid currently only supports the creation of structured grids with dimension == dimensionworld. Consider reading in a grid file instead.");
    }

    /*!
     * \brief Makes a structured cube grid using the structured grid factory
     */
    template <int dimension, int dimensionworld, std::enable_if_t<dimension == dimensionworld, int> = 0>
    void makeStructuredGrid(typename ParentType::CellType cellType,
                            const std::string& modelParamGroup)
    {
        // make a structured grid
        if (elType == Dune::cube)
            ParentType::template makeStructuredGrid<dimension, dimensionworld>(ParentType::CellType::Cube, modelParamGroup);
        else if (elType == Dune::simplex)
            ParentType::template makeStructuredGrid<dimension, dimensionworld>(ParentType::CellType::Simplex, modelParamGroup);
        else
            DUNE_THROW(Dune::IOError, "ALUGrid only supports Dune::cube or Dune::simplex as cell type!");
    }
};

/*!
 * \ingroup InputOutput
 * \brief Boundary flag
 */
//! alu uses boundary id
template<int dim, int dimworld, Dune::ALUGridElementType elType, Dune::ALUGridRefinementType refinementType>
class BoundaryFlag<Dune::ALUGrid<dim, dimworld, elType, refinementType>>
{
public:
    BoundaryFlag() : flag_(-1) {}

    template<class Intersection>
    BoundaryFlag(const Intersection& i) : flag_(-1)
    {
        if (i.boundary())
            flag_ = i.impl().boundaryId();
    }

    using value_type = int;

    value_type get() const { return flag_; }

private:
    int flag_;
};

#endif // HAVE_DUNE_ALUGRID

} // end namespace Dumux

#endif
