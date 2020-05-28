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
 * \ingroup BoxDFMModel
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format.
 */

#ifndef POROUSMEDIUMFLOW_BOXDFM_VTK_OUTPUT_MODULE_HH
#define POROUSMEDIUMFLOW_BOXDFM_VTK_OUTPUT_MODULE_HH

#include <set>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/io/vtkoutputmodule.hh>

namespace Dumux {

/*!
 * \ingroup BoxDFMModel
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format.
 *
 * This output module is specialized for writing out data obtained by the box-dfm
 * scheme. It writes out separate vtk files for the solution on the fracture. For
 * this, a grid type to be used the fracture-conforming lower-dimensional grid has
 * to be provided.
 *
 * \tparam TypeTag The TypeTag of the problem implementation
 * \tparam FractureGrid The Type used for the lower-dimensional grid
 *
 * Handles the output of scalar and vector fields to VTK formatted file for multiple
 * variables and time steps. Certain predefined fields can be registered on
 * initialization and/or be turned on/off using the designated properties. Additionally
 * non-standardized scalar and vector fields can be added to the writer manually.
 */
template<class GridVariables, class SolutionVector, class FractureGrid>
class BoxDfmVtkOutputModule : public VtkOutputModule<GridVariables, SolutionVector>
{
    using ParentType = VtkOutputModule<GridVariables, SolutionVector>;
    using GridGeometry = typename GridVariables::GridGeometry;
    using VV = typename GridVariables::VolumeVariables;
    using FluidSystem = typename VV::FluidSystem;
    using Scalar = typename GridVariables::Scalar;

    using GridView = typename GridGeometry::GridView;
    using FractureGridView = typename FractureGrid::LeafGridView;
    using FractureMapper = Dune::MultipleCodimMultipleGeomTypeMapper<FractureGridView>;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using GridIndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Field = Vtk::template Field<GridView>;
    using FractureField = Vtk::template Field<FractureGridView>;

    static_assert(dim > 1, "Box-Dfm output only works for dim > 1");
    static_assert(FractureGrid::dimension == int(dim-1), "Fracture grid must be of codimension one!");
    static_assert(FractureGrid::dimensionworld == int(dimWorld), "Fracture grid has to has the same coordinate dimension!");
    static_assert(GridGeometry::discMethod == DiscretizationMethod::box, "Box-Dfm output module can only be used with the box scheme!");
public:

    //! The constructor
    template< class FractureGridAdapter >
    BoxDfmVtkOutputModule(const GridVariables& gridVariables,
                          const SolutionVector& sol,
                          const std::string& name,
                          const FractureGridAdapter& fractureGridAdapter,
                          const std::string& paramGroup = "",
                          Dune::VTK::DataMode dm = Dune::VTK::conforming,
                          bool verbose = true)
    : ParentType(gridVariables, sol, name, paramGroup, dm, verbose)
    {
        // create the fracture grid and all objects needed on it
        initializeFracture_(fractureGridAdapter);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    //! Writing data
    //////////////////////////////////////////////////////////////////////////////////////////////

    //! Write the data for this timestep to file in four steps
    //! (1) We assemble all registered variable fields
    //! (2) We register them with the vtk writer
    //! (3) The writer writes the output for us
    //! (4) Clear the writer for the next time step
    void write(double time, Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
        Dune::Timer timer;

        // write to file depending on data mode
        const auto dm = this->dataMode();
        if (dm == Dune::VTK::conforming)
            writeConforming_(time, type);
        else if (dm == Dune::VTK::nonconforming)
            writeNonConforming_(time, type);
        else
            DUNE_THROW(Dune::NotImplemented, "Output for provided vtk data mode");

        //! output
        timer.stop();
        if (this->verbose())
            std::cout << "Writing output for problem \"" << this->name() << "\". Took " << timer.elapsed() << " seconds." << std::endl;
    }

private:
    //! Assembles the fields and adds them to the writer (conforming output)
    void writeConforming_(double time, Dune::VTK::OutputType type)
    {
        //////////////////////////////////////////////////////////////
        //! (1) Assemble all variable fields and add to writer
        //////////////////////////////////////////////////////////////

        // instatiate the velocity output
        std::vector<typename ParentType::VelocityOutput::VelocityVector> velocity;

        // process rank
        static bool addProcessRank = getParamFromGroup<bool>(this->paramGroup(), "Vtk.AddProcessRank");
        std::vector<double> rank;

        // volume variable data
        std::vector<std::vector<Scalar>> volVarScalarData;
        std::vector<std::vector<Scalar>> volVarScalarDataFracture;
        std::vector<std::vector<GlobalPosition>> volVarVectorData;
        std::vector<std::vector<GlobalPosition>> volVarVectorDataFracture;

        // some references for convenience
        const auto& gridView = this->gridGeometry().gridView();
        const auto& fractureGridView = fractureGrid_->leafGridView();
        const auto& volVarScalarDataInfo = this->volVarScalarDataInfo();
        const auto& volVarVectorDataInfo = this->volVarVectorDataInfo();

        //! Abort if no data was registered
        if (!volVarScalarDataInfo.empty()
            || !volVarVectorDataInfo.empty()
            || !this->fields().empty()
            || this->velocityOutput().enableOutput()
            || addProcessRank)
        {
            const auto numCells = gridView.size(0);
            const auto numDofs = gridView.size(dim);
            const auto numFractureVert = fractureGridView.size(FractureGridView::dimension);

            // get fields for all volume variables
            if (!this->volVarScalarDataInfo().empty())
            {
                volVarScalarData.resize(volVarScalarDataInfo.size(), std::vector<Scalar>(numDofs));
                volVarScalarDataFracture.resize(volVarScalarDataInfo.size(), std::vector<Scalar>(numFractureVert));
            }
            if (!this->volVarVectorDataInfo().empty())
            {
                volVarVectorData.resize(volVarVectorDataInfo.size(), std::vector<GlobalPosition>(numDofs));
                volVarVectorDataFracture.resize(volVarVectorDataInfo.size(), std::vector<GlobalPosition>(numFractureVert));
            }

            if (this->velocityOutput().enableOutput())
                for (int phaseIdx = 0; phaseIdx < this->velocityOutput().numFluidPhases(); ++phaseIdx)
                    velocity[phaseIdx].resize(numDofs);

            // maybe allocate space for the process rank
            if (addProcessRank) rank.resize(numCells);

            for (const auto& element : elements(gridView, Dune::Partitions::interior))
            {
                const auto eIdxGlobal = this->gridGeometry().elementMapper().index(element);

                auto fvGeometry = localView(this->gridGeometry());
                auto elemVolVars = localView(this->gridVariables().curGridVolVars());

                // If velocity output is enabled we need to bind to the whole stencil
                // otherwise element-local data is sufficient
                if (this->velocityOutput().enableOutput())
                {
                    fvGeometry.bind(element);
                    elemVolVars.bind(element, fvGeometry, this->sol());
                }
                else
                {
                    fvGeometry.bindElement(element);
                    elemVolVars.bindElement(element, fvGeometry, this->sol());
                }

                if (!volVarScalarDataInfo.empty() || !volVarVectorDataInfo.empty())
                {
                    for (auto&& scv : scvs(fvGeometry))
                    {
                        const auto dofIdxGlobal = scv.dofIndex();
                        const auto& volVars = elemVolVars[scv];

                        if (!scv.isOnFracture())
                        {
                            for (std::size_t i = 0; i < volVarScalarDataInfo.size(); ++i)
                                volVarScalarData[i][dofIdxGlobal] = volVarScalarDataInfo[i].get(volVars);
                            for (std::size_t i = 0; i < volVarVectorDataInfo.size(); ++i)
                                volVarVectorData[i][dofIdxGlobal] = volVarVectorDataInfo[i].get(volVars);
                        }
                        else
                        {
                            for (std::size_t i = 0; i < volVarScalarDataInfo.size(); ++i)
                                volVarScalarDataFracture[i][vertexToFractureVertexIdx_[dofIdxGlobal]] = volVarScalarDataInfo[i].get(volVars);
                            for (std::size_t i = 0; i < volVarVectorDataInfo.size(); ++i)
                                volVarVectorDataFracture[i][vertexToFractureVertexIdx_[dofIdxGlobal]] = volVarVectorDataInfo[i].get(volVars);
                        }
                    }
                }

                // velocity output
                if (this->velocityOutput().enableOutput())
                {
                    auto elemFluxVarsCache = localView(this->gridVariables().gridFluxVarsCache());
                    elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

                    for (int phaseIdx = 0; phaseIdx < this->velocityOutput().numFluidPhases(); ++phaseIdx)
                        this->velocityOutput().calculateVelocity(velocity[phaseIdx], element, fvGeometry, elemVolVars, elemFluxVarsCache, phaseIdx);
                }

                //! the rank
                if (addProcessRank)
                    rank[eIdxGlobal] = static_cast<double>(gridView.comm().rank());
            }

            //////////////////////////////////////////////////////////////
            //! (2) Register data fields with the vtk writer
            //////////////////////////////////////////////////////////////

            // volume variables if any
            for (std::size_t i = 0; i < volVarScalarDataInfo.size(); ++i)
            {
                this->sequenceWriter().addVertexData(volVarScalarData[i], volVarScalarDataInfo[i].name);
                fractureSequenceWriter_->addVertexData(volVarScalarDataFracture[i], volVarScalarDataInfo[i].name);
            }

            for (std::size_t i = 0; i < volVarVectorDataInfo.size(); ++i)
            {
                this->sequenceWriter().addVertexData( Field(gridView, this->gridGeometry().vertexMapper(), volVarVectorData[i],
                                                            volVarVectorDataInfo[i].name, /*numComp*/dimWorld, /*codim*/dim).get() );
                fractureSequenceWriter_->addVertexData( FractureField(fractureGridView, *fractureVertexMapper_, volVarVectorDataFracture[i],
                                                                      volVarVectorDataInfo[i].name, /*numComp*/dimWorld, /*codim*/dim-1).get() );
            }

            // the velocity field
            if (this->velocityOutput().enableOutput())
            {
                for (int phaseIdx = 0; phaseIdx < this->velocityOutput().numFluidPhases(); ++phaseIdx)
                    this->sequenceWriter().addVertexData( Field(gridView, this->gridGeometry().vertexMapper(), velocity[phaseIdx],
                                                                "velocity_" + std::string(this->velocityOutput().phaseName(phaseIdx)) + " (m/s)",
                                                                /*numComp*/dimWorld, /*codim*/dim).get() );
            }

            // the process rank
            if (addProcessRank)
                this->sequenceWriter().addCellData( Field(gridView, this->gridGeometry().elementMapper(), rank,
                                                          "process rank", /*numComp*/1, /*codim*/0).get() );

            // also register additional (non-standardized) user fields if any (only on matrix grid)
            for (auto&& field : this->fields())
            {
                if (field.codim() == 0)
                    this->sequenceWriter().addCellData(field.get());
                else if (field.codim() == dim)
                    this->sequenceWriter().addVertexData(field.get());
                else
                    DUNE_THROW(Dune::RangeError, "Cannot add wrongly sized vtk scalar field!");
            }
        }

        //////////////////////////////////////////////////////////////
        //! (2) The writer writes the output for us
        //////////////////////////////////////////////////////////////
        this->sequenceWriter().write(time, type);
        fractureSequenceWriter_->write(time, type);

        //////////////////////////////////////////////////////////////
        //! (3) Clear the writer
        //////////////////////////////////////////////////////////////
        this->writer().clear();
        fractureWriter_->clear();
    }

    //! Assembles the fields and adds them to the writer (conforming output)
    void writeNonConforming_(double time, Dune::VTK::OutputType type)
    {
        //////////////////////////////////////////////////////////////
        //! (1) Assemble all variable fields and add to writer
        //////////////////////////////////////////////////////////////

        // instatiate the velocity output
        std::vector<typename ParentType::VelocityOutput::VelocityVector> velocity;

        // process rank
        static bool addProcessRank = getParamFromGroup<bool>(this->paramGroup(), "Vtk.AddProcessRank");
        std::vector<double> rank;

        // volume variable data (indexing: volvardata/element/localcorner)
        using ScalarDataContainer = std::vector< std::vector<Scalar> >;
        using VectorDataContainer = std::vector< std::vector<GlobalPosition> >;
        std::vector< ScalarDataContainer > volVarScalarData;
        std::vector< ScalarDataContainer > volVarScalarDataFracture;
        std::vector< VectorDataContainer > volVarVectorData;
        std::vector< VectorDataContainer > volVarVectorDataFracture;

        // some references for convenience
        const auto& gridView = this->gridGeometry().gridView();
        const auto& fractureGridView = fractureGrid_->leafGridView();
        const auto& volVarScalarDataInfo = this->volVarScalarDataInfo();
        const auto& volVarVectorDataInfo = this->volVarVectorDataInfo();

        //! Abort if no data was registered
        if (!volVarScalarDataInfo.empty()
            || !volVarVectorDataInfo.empty()
            || !this->fields().empty()
            || this->velocityOutput().enableOutput()
            || addProcessRank)
        {
            const auto numCells = gridView.size(0);
            const auto numDofs = gridView.size(dim);
            const auto numFractureCells = fractureGridView.size(0);

            // get fields for all volume variables
            if (!this->volVarScalarDataInfo().empty())
            {
                volVarScalarData.resize(volVarScalarDataInfo.size(), ScalarDataContainer(numCells));
                volVarScalarDataFracture.resize(volVarScalarDataInfo.size(), ScalarDataContainer(numFractureCells));
            }
            if (!this->volVarVectorDataInfo().empty())
            {
                volVarVectorData.resize(volVarVectorDataInfo.size(), VectorDataContainer(numCells));
                volVarVectorDataFracture.resize(volVarVectorDataInfo.size(), VectorDataContainer(numFractureCells));
            }

            if (this->velocityOutput().enableOutput())
                for (int phaseIdx = 0; phaseIdx < this->velocityOutput().numFluidPhases(); ++phaseIdx)
                    velocity[phaseIdx].resize(numDofs);

            // maybe allocate space for the process rank
            if (addProcessRank) rank.resize(numCells);

            for (const auto& element : elements(gridView, Dune::Partitions::interior))
            {
                const auto eIdxGlobal = this->gridGeometry().elementMapper().index(element);
                const auto numCorners = element.subEntities(dim);

                auto fvGeometry = localView(this->gridGeometry());
                auto elemVolVars = localView(this->gridVariables().curGridVolVars());

                // resize element-local data containers (for bulk grid)
                for (std::size_t i = 0; i < volVarScalarDataInfo.size(); ++i)
                    volVarScalarData[i][eIdxGlobal].resize(numCorners);
                for (std::size_t i = 0; i < volVarVectorDataInfo.size(); ++i)
                    volVarVectorData[i][eIdxGlobal].resize(numCorners);

                // If velocity output is enabled we need to bind to the whole stencil
                // otherwise element-local data is sufficient
                if (this->velocityOutput().enableOutput())
                {
                    fvGeometry.bind(element);
                    elemVolVars.bind(element, fvGeometry, this->sol());
                }
                else
                {
                    fvGeometry.bindElement(element);
                    elemVolVars.bindElement(element, fvGeometry, this->sol());
                }

                if (!volVarScalarDataInfo.empty() || !volVarVectorDataInfo.empty())
                {
                    for (auto&& scv : scvs(fvGeometry))
                    {
                        const auto& volVars = elemVolVars[scv];

                        if (!scv.isOnFracture())
                        {
                            for (std::size_t i = 0; i < volVarScalarDataInfo.size(); ++i)
                                volVarScalarData[i][eIdxGlobal][scv.localDofIndex()] = volVarScalarDataInfo[i].get(volVars);
                            for (std::size_t i = 0; i < volVarVectorDataInfo.size(); ++i)
                                volVarVectorData[i][eIdxGlobal][scv.localDofIndex()]  = volVarVectorDataInfo[i].get(volVars);
                        }
                        else
                        {
                            const auto fIdx = scv.facetIndexInElement();
                            const auto& localMap = fractureElementMap_[eIdxGlobal];
                            const auto fracEIdx = std::find_if(localMap.begin(), localMap.end(), [fIdx] (const auto& p) { return p.first == fIdx; })->second;
                            for (std::size_t i = 0; i < volVarScalarDataInfo.size(); ++i)
                                volVarScalarDataFracture[i][fracEIdx].push_back(volVarScalarDataInfo[i].get(volVars));
                            for (std::size_t i = 0; i < volVarVectorDataInfo.size(); ++i)
                                volVarVectorDataFracture[i][fracEIdx].push_back(volVarVectorDataInfo[i].get(volVars));
                        }
                    }
                }

                // velocity output
                if (this->velocityOutput().enableOutput())
                {
                    auto elemFluxVarsCache = localView(this->gridVariables().gridFluxVarsCache());
                    elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

                    for (int phaseIdx = 0; phaseIdx < this->velocityOutput().numFluidPhases(); ++phaseIdx)
                        this->velocityOutput().calculateVelocity(velocity[phaseIdx], element, fvGeometry, elemVolVars, elemFluxVarsCache, phaseIdx);
                }

                //! the rank
                if (addProcessRank)
                    rank[eIdxGlobal] = static_cast<double>(gridView.comm().rank());
            }

            //////////////////////////////////////////////////////////////
            //! (2) Register data fields with the vtk writer
            //////////////////////////////////////////////////////////////

            // volume variables if any
            for (std::size_t i = 0; i < volVarScalarDataInfo.size(); ++i)
            {
                this->sequenceWriter().addVertexData( Field(gridView, this->gridGeometry().elementMapper(), volVarScalarData[i],
                                                            volVarScalarDataInfo[i].name, /*numComp*/1, /*codim*/dim,
                                                            /*nonconforming*/this->dataMode()).get() );
                fractureSequenceWriter_->addVertexData( FractureField(fractureGridView, *fractureElementMapper_, volVarScalarDataFracture[i],
                                                                      volVarScalarDataInfo[i].name, /*numComp*/1, /*codim*/dim-1,
                                                                      /*nonconforming*/this->dataMode()).get() );
            }

            for (std::size_t i = 0; i < volVarVectorDataInfo.size(); ++i)
            {
                this->sequenceWriter().addVertexData( Field(gridView, this->gridGeometry().elementMapper(), volVarVectorData[i],
                                                            volVarVectorDataInfo[i].name, /*numComp*/dimWorld, /*codim*/dim,
                                                            /*nonconforming*/this->dataMode()).get() );
                fractureSequenceWriter_->addVertexData( FractureField(fractureGridView, *fractureElementMapper_, volVarVectorDataFracture[i],
                                                                      volVarVectorDataInfo[i].name, /*numComp*/dimWorld, /*codim*/dim-1,
                                                                      /*nonconforming*/this->dataMode()).get() );
            }

            // the velocity field
            if (this->velocityOutput().enableOutput())
            {
                for (int phaseIdx = 0; phaseIdx < this->velocityOutput().numFluidPhases(); ++phaseIdx)
                    this->sequenceWriter().addVertexData( Field(gridView, this->gridGeometry().vertexMapper(), velocity[phaseIdx],
                                                                "velocity_" + std::string(this->velocityOutput().phaseName(phaseIdx)) + " (m/s)",
                                                                /*numComp*/dimWorld, /*codim*/dim).get() );
            }

            // the process rank
            if (addProcessRank)
                this->sequenceWriter().addCellData( Field(gridView, this->gridGeometry().elementMapper(), rank,
                                                          "process rank", /*numComp*/1, /*codim*/0).get() );

            // also register additional (non-standardized) user fields if any (only on matrix grid)
            for (auto&& field : this->fields())
            {
                if (field.codim() == 0)
                    this->sequenceWriter().addCellData(field.get());
                else if (field.codim() == dim)
                    this->sequenceWriter().addVertexData(field.get());
                else
                    DUNE_THROW(Dune::RangeError, "Cannot add wrongly sized vtk scalar field!");
            }
        }

        //////////////////////////////////////////////////////////////
        //! (2) The writer writes the output for us
        //////////////////////////////////////////////////////////////
        this->sequenceWriter().write(time, type);
        fractureSequenceWriter_->write(time, type);

        //////////////////////////////////////////////////////////////
        //! (3) Clear the writer
        //////////////////////////////////////////////////////////////
        this->writer().clear();
        fractureWriter_->clear();
    }

    //! Creates the lower-dimensional fracture grid, index maps and writers
    template< class FractureGridAdapter >
    void initializeFracture_(const FractureGridAdapter& fractureGridAdapter)
    {
        const auto& gridGeometry = this->gridGeometry();
        const auto& gridView = gridGeometry.gridView();
        Dune::GridFactory<FractureGrid> gridFactory;

        // insert fracture vertices
        std::size_t fracVertexCount = 0;
        vertexToFractureVertexIdx_.resize(gridView.size(dim));
        for (const auto& v : vertices(gridView))
        {
            if (fractureGridAdapter.isOnFacetGrid(v))
            {
                gridFactory.insertVertex(v.geometry().center());
                vertexToFractureVertexIdx_[gridGeometry.vertexMapper().index(v)] = fracVertexCount++;
            }
        }

        // insert fracture elements
        std::size_t fractureElementCount = 0;
        fractureElementMap_.resize(gridView.size(0));
        std::set< std::pair<GridIndexType, unsigned int> > handledFacets;
        for (const auto& element : elements(gridView))
        {
            const auto eIdxGlobal = gridGeometry.elementMapper().index(element);
            const auto refElement = referenceElement(element);

            for (const auto& is : intersections(gridView, element))
            {
                // obtain all vertex indices on this intersection
                const auto& isGeometry = is.geometry();
                const auto numCorners = isGeometry.corners();
                const auto indexInInside = is.indexInInside();

                std::vector<GridIndexType> isVertexIndices(numCorners);
                for (unsigned int i = 0; i < numCorners; ++i)
                    isVertexIndices[i] = gridGeometry.vertexMapper().subIndex(element,
                                                                                refElement.subEntity(indexInInside, 1, i, dim),
                                                                                dim);

                // determine if this is a fracture facet & if it has to be inserted
                bool insertFacet = false;
                if (fractureGridAdapter.composeFacetElement(isVertexIndices))
                {
                    insertFacet = true;
                    if (!is.boundary())
                    {
                        // only proceed if facet has not been handled yet
                        const auto outsideEIdx = gridGeometry.elementMapper().index(is.outside());
                        const auto idxInOutside = is.indexInOutside();
                        const auto pair = std::make_pair(outsideEIdx, idxInOutside);
                        if (handledFacets.count( pair ) != 0)
                        {
                            insertFacet = false;

                            // obtain the fracture grid elem idx from outside map and insert to map
                            const auto& outsideMap = fractureElementMap_[outsideEIdx];
                            auto it = std::find_if(outsideMap.begin(), outsideMap.end(), [idxInOutside] (const auto& p) { return p.first == idxInOutside; });
                            fractureElementMap_[eIdxGlobal].push_back( std::make_pair(indexInInside, it->second) );
                        }
                    }
                }

                if (insertFacet)
                {
                    // transform intersection vertex indices to frac grid indices
                    std::for_each( isVertexIndices.begin(),
                                   isVertexIndices.end(),
                                   [&] (auto& idx) { idx = this->vertexToFractureVertexIdx_[idx]; } );

                    // insert the element
                    gridFactory.insertElement(isGeometry.type(), isVertexIndices);

                    // insert to set of handled facets
                    handledFacets.insert( std::make_pair(eIdxGlobal, indexInInside) );
                    fractureElementMap_[eIdxGlobal].push_back( std::make_pair(indexInInside, fractureElementCount) );
                    fractureElementCount++;
                }
            }
        }

        // make grid and get grid view
        fractureGrid_ = std::shared_ptr<FractureGrid>(gridFactory.createGrid());

        // update fracture mappers
        const auto& fractureGridView = fractureGrid_->leafGridView();
        fractureVertexMapper_ = std::make_unique<FractureMapper>(fractureGridView, Dune::mcmgVertexLayout());
        fractureElementMapper_ = std::make_unique<FractureMapper>(fractureGridView, Dune::mcmgElementLayout());

        // obtain map fracture insertion indices -> fracture grid indices
        std::vector<GridIndexType> insToVertexIdx(fractureGridView.size(FractureGridView::dimension));
        std::vector<GridIndexType> insToElemIdx(fractureGridView.size(0));
        for (const auto& v : vertices(fractureGridView)) insToVertexIdx[ gridFactory.insertionIndex(v) ] = fractureVertexMapper_->index(v);
        for (const auto& e : elements(fractureGridView)) insToElemIdx[ gridFactory.insertionIndex(e) ] = fractureElementMapper_->index(e);

        // update vertex index map
        for (GridIndexType dofIdx = 0; dofIdx < gridView.size(GridView::dimension); ++dofIdx)
            if (gridGeometry.dofOnFracture(dofIdx))
                vertexToFractureVertexIdx_[dofIdx] = insToVertexIdx[ vertexToFractureVertexIdx_[dofIdx] ];

        // update fracture element map
        for (auto& elemLocalMap : fractureElementMap_)
            for (auto& dataPair : elemLocalMap)
                dataPair.second = insToElemIdx[ dataPair.second ];

        // instantiate writers for the fracture
        fractureWriter_ = std::make_shared< Dune::VTKWriter<FractureGridView> >(fractureGridView, this->dataMode());
        fractureSequenceWriter_ = std::make_unique< Dune::VTKSequenceWriter<FractureGridView> >(fractureWriter_, this->name() + "_fracture");
    }

    std::shared_ptr<FractureGrid> fractureGrid_;

    std::unique_ptr<FractureMapper> fractureVertexMapper_;
    std::unique_ptr<FractureMapper> fractureElementMapper_;

    std::shared_ptr<Dune::VTKWriter<FractureGridView>> fractureWriter_;
    std::unique_ptr< Dune::VTKSequenceWriter<FractureGridView> > fractureSequenceWriter_;

    // maps to a bulk grid vertex the vertex index within the fracture grid
    std::vector<GridIndexType> vertexToFractureVertexIdx_;

    // maps to the local facet indices of an element the corresponding fracture element indices
    std::vector< std::vector<std::pair<GridIndexType, unsigned int>> > fractureElementMap_;
};

} // end namespace Dumux

#endif
