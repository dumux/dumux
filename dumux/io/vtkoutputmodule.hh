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
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 */
#ifndef VTK_OUTPUT_MODULE_HH
#define VTK_OUTPUT_MODULE_HH

#include <dune/common/fvector.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include <dumux/io/vtknestedfunction.hh>

namespace Properties
{
NEW_PROP_TAG(VtkAddVelocity);
NEW_PROP_TAG(VtkAddProcessRank);
NEW_PROP_TAG(VtkAddPorosity);
NEW_PROP_TAG(VtkAddPermeability);
}

namespace Dumux
{

/*!
 * \ingroup InputOutput
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 *
 * Handles the output of scalar and vector fields to VTK formatted file for multiple
 * variables and timesteps. Certain predefined fields can be registered on problem / model
 * initialization and/or be turned on/off using the designated properties. Additionally
 * non-standardized scalar and vector fields can be added to the writer manually.
 */
template<typename TypeTag>
class VtkOutputModule
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ElementMapper = typename GET_PROP_TYPE(TypeTag, ElementMapper);
    using VertexMapper = typename GET_PROP_TYPE(TypeTag, VertexMapper);

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);

    struct PriVarScalarDataInfo { unsigned int pvIdx; std::string name; };
    struct PriVarVectorDataInfo { std::vector<unsigned int> pvIdx; std::string name; };
    struct SecondVarScalarDataInfo { std::function<Scalar(const VolumeVariables&)> get; std::string name; };

public:

    VtkOutputModule(const Problem& problem,
                    Dune::VTK::DataMode dm = Dune::VTK::conforming)
    : problem_(problem),
      writer_(std::make_shared<Dune::VTKWriter<GridView>>(problem.gridView(), dm)),
      sequenceWriter_(writer_, problem.name())
    {}

    //////////////////////////////////////////////////////////////////////////////////////////////
    //! Methods to conveniently add primary and secondary variables upon problem initialization
    //! Do not call these methods after initialization
    //////////////////////////////////////////////////////////////////////////////////////////////

    //! Output a scalar primary variable
    //! \param name The name of the vtk field
    //! \param pvIdx The index in the primary variables vector
    void addPrimaryVariable(const std::string& name, unsigned int pvIdx)
    {
        priVarScalarDataInfo_.push_back(PriVarScalarDataInfo{pvIdx, name});
    }

    //! Output a vector primary variable
    //! \param name The name of the vtk field
    //! \param pvIndices A vector of indices in the primary variables vector to group for vector visualization
    void addPrimaryVariable(const std::string& name, std::vector<unsigned int> pvIndices)
    {
        assert(pvIndices.size() < 4 && "Vtk doesn't support vector dimensions greater than 3!");
        priVarVectorDataInfo_.push_back(PriVarVectorDataInfo{pvIndices, name});
    }

    //! Output a secondary scalar variable
    //! \param name The name of the vtk field
    //! \param f A function taking a VolumeVariables object and returning the desired scalar
    void addSecondaryVariable(const std::string& name, std::function<Scalar(const VolumeVariables&)>&& f)
    {
        secondVarScalarDataInfo_.push_back(SecondVarScalarDataInfo{f, name});
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    //! Methods to add addtional (non-standardized) scalar or vector fields to the output queue
    //! Call these methods on the output module obtained as argument of the addVtkOutputFields
    //! method that can be overloaded in the problem implementation
    //////////////////////////////////////////////////////////////////////////////////////////////

    //! Output a scalar field
    //! \param name The name of the vtk field
    //! \param codim The codimension of the entities the vtk field lives on (options: 0 := elements, dim := vertices)
    //! \returns A reference to the resized scalar field to be filled with the actual data
    std::vector<Scalar>& createScalarField(const std::string& name, int codim)
    {
        scalarFields_.emplace_back(std::make_pair(std::vector<Scalar>(problem_.gridView().size(codim)), name));
        return scalarFields_.back().first;
    }

    //! Output a vector field
    //! \param name The name of the vtk field
    //! \param codim The codimension of the entities the vtk field lives on (options: 0 := elements, dim := vertices)
    //! \returns A reference to the resized vector field to be filled with the actual data
    std::vector<GlobalPosition>& createVectorField(const std::string& name, int codim)
    {
        vectorFields_.emplace_back(std::make_pair(std::vector<GlobalPosition>(problem_.gridView().size(codim)), name));
        return vectorFields_.back().first;
    }

    //! Write the data for this timestep to file in four steps
    //! (1) Allow user to register additional (non-standardized) fields
    //! (2) We assemble all registered variable fields
    //! (3) We register them with the vtk writer
    //! (4) The writer writes the output for us
    //! (5) Clear the writer for the next time step
    void write(double time, Dune::VTK::OutputType type = Dune::VTK::ascii)
    {

        //////////////////////////////////////////////////////////////
        //! (1) Register addtional (non-standardized) data fields with the vtk writer
        //!     Using the add scalar field or vector field methods
        //////////////////////////////////////////////////////////////
        problem_.model().addVtkOutputFields(*this);
        problem_.addVtkOutputFields(*this);

        //! Abort if no data was registered
        //! \todo This is not necessary anymore once the old style multiwriter is removed
        if (priVarScalarDataInfo_.empty()
            && priVarVectorDataInfo_.empty()
            && secondVarScalarDataInfo_.empty()
            && scalarFields_.empty()
            && vectorFields_.empty())
            return;

        //////////////////////////////////////////////////////////////
        //! (2) Assemble all variable fields with registered info
        //////////////////////////////////////////////////////////////
        auto numCells = problem_.gridView().size(0);
        auto numDofs = problem_.model().numDofs();

        // get fields for all primary variables
        std::vector<std::vector<Scalar>> priVarScalarData(priVarScalarDataInfo_.size());
        for (auto&& p : priVarScalarData)
            p.resize(numDofs);

        std::vector<std::vector<Scalar>> priVarVectorData(priVarVectorDataInfo_.size());
        for (std::size_t i = 0; i < priVarVectorDataInfo_.size(); ++i)
            priVarVectorData[i].resize(numDofs*priVarVectorDataInfo_[i].pvIdx.size());

        // get fields for all secondary variables
        std::vector<std::vector<Scalar>> secondVarScalarData(secondVarScalarDataInfo_.size());
        for (auto&& s : secondVarScalarData)
            s.resize(numDofs);

        // instatiate the velocity output
        ImplicitVelocityOutput<TypeTag> velocityOutput(problem_);
        std::array<std::vector<GlobalPosition>, numPhases> velocity;

        if (velocityOutput.enableOutput())
        {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                velocity[phaseIdx].resize(numDofs);
        }

        // maybe allocate space for the process rank
        std::vector<Scalar> rank;
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddProcessRank))
            rank.resize(numCells);

        std::vector<Scalar> porosity;
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddPorosity))
            porosity.resize(numDofs);

        std::vector<Scalar> permeability;
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddPermeability))
            permeability.resize(numDofs);

        for (const auto& element : elements(problem_.gridView(), Dune::Partitions::interior))
        {
            const auto eIdxGlobal = problem_.elementMapper().index(element);

            // cell-centered models
            if(!isBox)
            {
                //! primary variable data
                for (std::size_t i = 0; i < priVarScalarDataInfo_.size(); ++i)
                    priVarScalarData[i][eIdxGlobal] = problem_.model().curSol()[eIdxGlobal][priVarScalarDataInfo_[i].pvIdx];

                for (std::size_t i = 0; i < priVarVectorDataInfo_.size(); ++i)
                    for (std::size_t j = 0; j < priVarVectorDataInfo_[i].pvIdx.size(); ++j)
                        priVarVectorData[i][eIdxGlobal*priVarVectorDataInfo_[i].pvIdx.size() + j]
                            = problem_.model().curSol()[eIdxGlobal][priVarVectorDataInfo_[i].pvIdx[j]];
            }

            auto fvGeometry = localView(problem_.model().globalFvGeometry());
            auto elemVolVars = localView(problem_.model().curGlobalVolVars());
            auto curElemSol = problem_.model().elementSolution(element, problem_.model().curSol());

            // If velocity output is enabled we need to bind to the whole stencil
            // otherwise element-local data is sufficient
            if (velocityOutput.enableOutput())
                fvGeometry.bind(element);
            else
                fvGeometry.bindElement(element);

            // If velocity output is enabled we need to bind to the whole stencil
            // otherwise element-local data is sufficient
            if (velocityOutput.enableOutput())
                elemVolVars.bind(element, fvGeometry, problem_.model().curSol());
            else
                elemVolVars.bindElement(element, fvGeometry, problem_.model().curSol());

            for (auto&& scv : scvs(fvGeometry))
            {
                const auto dofIdxGlobal = scv.dofIndex();

                // for box model do the privars here
                if (isBox)
                {
                    //! primary variable data
                    for (std::size_t i = 0; i < priVarScalarDataInfo_.size(); ++i)
                        priVarScalarData[i][dofIdxGlobal] = problem_.model().curSol()[dofIdxGlobal][priVarScalarDataInfo_[i].pvIdx];

                    for (std::size_t i = 0; i < priVarVectorDataInfo_.size(); ++i)
                        for (std::size_t j = 0; j < priVarVectorDataInfo_[i].pvIdx.size(); ++j)
                            priVarVectorData[i][dofIdxGlobal*priVarVectorDataInfo_[i].pvIdx.size() + j]
                                = problem_.model().curSol()[dofIdxGlobal][priVarVectorDataInfo_[i].pvIdx[j]];

                }

                // secondary variables
                const auto& volVars = elemVolVars[scv];

                for (std::size_t i = 0; i < secondVarScalarDataInfo_.size(); ++i)
                    secondVarScalarData[i][dofIdxGlobal] = secondVarScalarDataInfo_[i].get(volVars);

                if (GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddPorosity))
                    porosity[dofIdxGlobal] = problem_.spatialParams().porosity(element, scv, curElemSol);

                if (GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddPermeability))
                    permeability[dofIdxGlobal] = problem_.spatialParams().permeability(element, scv, curElemSol);

            }

            // velocity output
            if (velocityOutput.enableOutput())
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    velocityOutput.calculateVelocity(velocity[phaseIdx], elemVolVars, fvGeometry, element, phaseIdx);

            //! the rank
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddProcessRank))
                rank[eIdxGlobal] = problem_.gridView().comm().rank();
        }

        //////////////////////////////////////////////////////////////
        //! (3) Register data fields with the vtk writer
        //////////////////////////////////////////////////////////////

        // primary variables
        for (std::size_t i = 0; i < priVarScalarDataInfo_.size(); ++i)
            addDofDataForWriter_(sequenceWriter_, priVarScalarData[i], priVarScalarDataInfo_[i].name);

        for (std::size_t i = 0; i < priVarVectorDataInfo_.size(); ++i)
            addDofDataForWriter_(sequenceWriter_, priVarVectorData[i], priVarVectorDataInfo_[i].name, priVarVectorDataInfo_[i].pvIdx.size());

        // secondary variables
        for (std::size_t i = 0; i < secondVarScalarDataInfo_.size(); ++i)
            addDofDataForWriter_(sequenceWriter_, secondVarScalarData[i], secondVarScalarDataInfo_[i].name);

        // the velocity field
        if (velocityOutput.enableOutput())
        {
            if (isBox)
            {
                using NestedFunction = VtkNestedFunction<GridView, VertexMapper, std::vector<GlobalPosition>>;
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    sequenceWriter_.addVertexData(std::make_shared<NestedFunction>("velocity_" + std::string(FluidSystem::phaseName(phaseIdx)),
                                                                                 problem_.gridView(), problem_.vertexMapper(),
                                                                                 velocity[phaseIdx], dim, dimWorld));
            }
            // cell-centered models
            else
            {
                using NestedFunction = VtkNestedFunction<GridView, ElementMapper, std::vector<GlobalPosition>>;
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    sequenceWriter_.addCellData(std::make_shared<NestedFunction>("velocity_" + std::string(FluidSystem::phaseName(phaseIdx)),
                                                                                 problem_.gridView(), problem_.elementMapper(),
                                                                                 velocity[phaseIdx], 0, dimWorld));
            }
        }

        // the process rank
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddProcessRank))
            sequenceWriter_.addCellData(rank, "process rank");

        // the porosity
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddPorosity))
            addDofDataForWriter_(sequenceWriter_, porosity, "porosity");

        // the permeability
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddPermeability))
            addDofDataForWriter_(sequenceWriter_, permeability, "permeability");

        // also register additional (non-standardized) user fields
        for (std::size_t i = 0; i < scalarFields_.size(); ++i)
        {
            if (scalarFields_[i].first.size() == std::size_t(problem_.gridView().size(0)))
                sequenceWriter_.addCellData(scalarFields_[i].first, scalarFields_[i].second);
            else if (scalarFields_[i].first.size() == std::size_t(problem_.gridView().size(dim)))
                sequenceWriter_.addVertexData(scalarFields_[i].first, scalarFields_[i].second);
            else
                DUNE_THROW(Dune::RangeError, "Cannot add wrongly sized vtk scalar field!");
        }

        for (std::size_t i = 0; i < vectorFields_.size(); ++i)
        {
            if (scalarFields_[i].first.size() == std::size_t(problem_.gridView().size(0)))
            {
                using NestedFunction = VtkNestedFunction<GridView, ElementMapper, std::vector<GlobalPosition>>;
                sequenceWriter_.addCellData(std::make_shared<NestedFunction>(vectorFields_[i].second,
                                                                             problem_.gridView(), problem_.elementMapper(),
                                                                             vectorFields_[i].first, 0, dimWorld));
            }
            else if (scalarFields_[i].first.size() == std::size_t(problem_.gridView().size(dim)))
            {
                using NestedFunction = VtkNestedFunction<GridView, VertexMapper, std::vector<GlobalPosition>>;
                sequenceWriter_.addVertexData(std::make_shared<NestedFunction>(vectorFields_[i].second,
                                                                               problem_.gridView(), problem_.vertexMapper(),
                                                                               vectorFields_[i].first, dim, dimWorld));
            }
            else
                DUNE_THROW(Dune::RangeError, "Cannot add wrongly sized vtk vector field!");
        }

        //////////////////////////////////////////////////////////////
        //! (4) The writer writes the output for us
        //////////////////////////////////////////////////////////////
        sequenceWriter_.write(time, type);

        //////////////////////////////////////////////////////////////
        //! (5) Clear all fields and writer for the next time step
        //////////////////////////////////////////////////////////////
        clear();
    }

    //! clear all data in the writer
    void clear()
    {
        writer_->clear();
        scalarFields_.clear();
        vectorFields_.clear();
    }

    /*!
     * \brief This method writes the complete state of the problem
     *        to the harddisk.
     *
     * The file will start with the prefix returned by the name()
     * method, has the current time of the simulation clock in it's
     * name and uses the extension <tt>.drs</tt>. (Dumux ReStart
     * file.)  See Restart for details.
     *
     * \tparam Restarter The serializer type
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        // TODO implement
    }

    /*!
     * \brief This method restores the complete state of the problem
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     *
     * \tparam Restarter The deserializer type
     *
     * \param res The deserializer object
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        // TODO implement
    }

private:

    template<typename Writer, typename... Args>
    void addDofDataForWriter_(Writer& writer,
                              Args&&... args)
    {
        if (isBox)
            writer.addVertexData(std::forward<Args>(args)...);
        else
            writer.addCellData(std::forward<Args>(args)...);
    }

    const Problem& problem_;
    std::shared_ptr<Dune::VTKWriter<GridView>> writer_;
    Dune::VTKSequenceWriter<GridView> sequenceWriter_;

    std::vector<PriVarScalarDataInfo> priVarScalarDataInfo_;
    std::vector<PriVarVectorDataInfo> priVarVectorDataInfo_;
    std::vector<SecondVarScalarDataInfo> secondVarScalarDataInfo_;

    std::vector<std::pair<std::vector<Scalar>, std::string>> scalarFields_;
    std::vector<std::pair<std::vector<GlobalPosition>, std::string>> vectorFields_;
};

} // end namespace Dumux

#endif
