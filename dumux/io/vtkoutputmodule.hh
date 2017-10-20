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

#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dumux/implicit/properties.hh>
#include <dumux/io/vtknestedfunction.hh>

namespace Properties
{
NEW_PROP_TAG(VtkAddVelocity);
NEW_PROP_TAG(VtkAddProcessRank);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(NumPhases);
}

namespace Dumux
{

namespace VtkImplDetail
{
    auto defaultAdderFunction = [](auto& m){};
}

/*!
 * \ingroup InputOutput
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 *
 * Handles the output of scalar and vector fields to VTK formatted file for multiple
 * variables and timesteps. Certain predefined fields can be registered on
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
    using ElementMapper = typename GET_PROP_TYPE(TypeTag, ElementMapper);
    using VertexMapper = typename GET_PROP_TYPE(TypeTag, VertexMapper);
    using Implementation = typename GET_PROP_TYPE(TypeTag, VtkOutputModule);
    using VelocityOutput = typename GET_PROP_TYPE(TypeTag, VelocityOutput);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);
    static constexpr int dofCodim = isBox ? dim : 0;

    struct PriVarScalarDataInfo { unsigned int pvIdx; std::string name; };
    struct PriVarVectorDataInfo { std::vector<unsigned int> pvIdx; std::string name; };
    struct SecondVarScalarDataInfo { std::function<Scalar(const VolumeVariables&)> get; std::string name; };

    const std::string modelParamGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);

public:

    VtkOutputModule(const Problem& problem,
                    const FVGridGeometry& fvGridGeometry,
                    const GridVariables& gridVariables,
                    const SolutionVector& sol,
                    const std::string& name,
                    bool verbose = true,
                    Dune::VTK::DataMode dm = Dune::VTK::conforming)
    : problem_(problem)
    , gridGeom_(fvGridGeometry)
    , gridVariables_(gridVariables)
    , sol_(sol)
    , name_(name)
    , verbose_(fvGridGeometry.gridView().comm().rank() == 0 && verbose)
    , writer_(std::make_shared<Dune::VTKWriter<GridView>>(fvGridGeometry.gridView(), dm))
    , sequenceWriter_(writer_, name)
    {}

    //////////////////////////////////////////////////////////////////////////////////////////////
    //! Methods to conveniently add primary and secondary variables upon initialization
    //! Do not call these methods after initialization i.e. not within the time loop
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
    //! Call these methods on the output module obtained as argument of the addVtkOutputFields function
    //////////////////////////////////////////////////////////////////////////////////////////////

    //! Output a scalar field
    //! \param name The name of the vtk field
    //! \param codim The codimension of the entities the vtk field lives on (options: 0 := elements, dim := vertices)
    //! \returns A reference to the resized scalar field to be filled with the actual data
    std::vector<Scalar>& createScalarField(const std::string& name, int codim)
    {
        scalarFields_.emplace_back(std::make_pair(std::vector<Scalar>(gridGeom_.gridView().size(codim)), name));
        return scalarFields_.back().first;
    }

    //! Output a vector field
    //! \param name The name of the vtk field
    //! \param codim The codimension of the entities the vtk field lives on (options: 0 := elements, dim := vertices)
    //! \returns A reference to the resized vector field to be filled with the actual data
    std::vector<GlobalPosition>& createVectorField(const std::string& name, int codim)
    {
        vectorFields_.emplace_back(std::make_pair(std::vector<GlobalPosition>(gridGeom_.gridView().size(codim)), name));
        return vectorFields_.back().first;
    }

    //! Write the data for this timestep to file in four steps
    //! (1) Allow user to register additional (non-standardized) fields
    //! (2) We assemble all registered variable fields
    //! (3) We register them with the vtk writer
    //! (4) The writer writes the output for us
    //! (5) Clear the writer for the next time step
    template<typename AdderFunction = decltype(VtkImplDetail::defaultAdderFunction)>
    void write(double time,
               const AdderFunction& addVtkOutputFields = VtkImplDetail::defaultAdderFunction,
               Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
        Dune::Timer timer;
        //////////////////////////////////////////////////////////////
        //! (1) Register addtional (non-standardized) data fields with the vtk writer
        //!     Using the add scalar field or vector field methods
        //////////////////////////////////////////////////////////////
        addVtkOutputFields(asImp_());

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
        auto numCells = gridGeom_.gridView().size(0);
        auto numDofs = asImp_().numDofs_();

        // get fields for all primary variables
        std::vector<std::vector<Scalar>> priVarScalarData(priVarScalarDataInfo_.size(), std::vector<Scalar>(numDofs));

        std::vector<std::vector<Scalar>> priVarVectorData(priVarVectorDataInfo_.size());
        for (std::size_t i = 0; i < priVarVectorDataInfo_.size(); ++i)
            priVarVectorData[i].resize(numDofs*priVarVectorDataInfo_[i].pvIdx.size());

        // get fields for all secondary variables
        std::vector<std::vector<Scalar>> secondVarScalarData(secondVarScalarDataInfo_.size(), std::vector<Scalar>(numDofs));

        // instatiate the velocity output
        VelocityOutput velocityOutput(problem_, gridGeom_, gridVariables_, sol_);
        std::array<std::vector<GlobalPosition>, numPhases> velocity;

        if (velocityOutput.enableOutput())
        {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                if(isBox && dim == 1)
                    velocity[phaseIdx].resize(numCells);
                else
                    velocity[phaseIdx].resize(numDofs);
            }
        }

        // maybe allocate space for the process rank
        std::vector<Scalar> rank;
        static bool addProcessRank = getParamFromGroup<bool>(modelParamGroup, "Vtk.AddProcessRank");
        if (addProcessRank) rank.resize(numCells);

        for (const auto& element : elements(gridGeom_.gridView(), Dune::Partitions::interior))
        {
            const auto eIdxGlobal = gridGeom_.elementMapper().index(element);

            // cell-centered models
            if(!isBox)
            {
                //! primary variable data
                for (std::size_t i = 0; i < priVarScalarDataInfo_.size(); ++i)
                    priVarScalarData[i][eIdxGlobal] = asImp_().getPriVarData_(eIdxGlobal, priVarScalarDataInfo_[i].pvIdx);

                for (std::size_t i = 0; i < priVarVectorDataInfo_.size(); ++i)
                    for (std::size_t j = 0; j < priVarVectorDataInfo_[i].pvIdx.size(); ++j)
                        priVarVectorData[i][eIdxGlobal*priVarVectorDataInfo_[i].pvIdx.size() + j]
                            = asImp_().getPriVarData_(eIdxGlobal, priVarVectorDataInfo_[i].pvIdx[j]);
            }

            auto fvGeometry = localView(gridGeom_);
            auto elemVolVars = localView(gridVariables_.curGridVolVars());

            // If velocity output is enabled we need to bind to the whole stencil
            // otherwise element-local data is sufficient
            if (velocityOutput.enableOutput())
                fvGeometry.bind(element);
            else
                fvGeometry.bindElement(element);

            // If velocity output is enabled we need to bind to the whole stencil
            // otherwise element-local data is sufficient
            if (velocityOutput.enableOutput())
                elemVolVars.bind(element, fvGeometry, sol_);
            else
                elemVolVars.bindElement(element, fvGeometry, sol_);

            for (auto&& scv : scvs(fvGeometry))
            {
                const auto dofIdxGlobal = scv.dofIndex();

                // for box model do the privars here
                if (isBox)
                {
                    //! primary variable data
                    for (std::size_t i = 0; i < priVarScalarDataInfo_.size(); ++i)
                        priVarScalarData[i][dofIdxGlobal] = asImp_().getPriVarData_(dofIdxGlobal, priVarScalarDataInfo_[i].pvIdx);

                    for (std::size_t i = 0; i < priVarVectorDataInfo_.size(); ++i)
                        for (std::size_t j = 0; j < priVarVectorDataInfo_[i].pvIdx.size(); ++j)
                            priVarVectorData[i][dofIdxGlobal*priVarVectorDataInfo_[i].pvIdx.size() + j]
                                = asImp_().getPriVarData_(dofIdxGlobal,priVarVectorDataInfo_[i].pvIdx[j]);

                }

                // secondary variables
                const auto& volVars = elemVolVars[scv];

                for (std::size_t i = 0; i < secondVarScalarDataInfo_.size(); ++i)
                    secondVarScalarData[i][dofIdxGlobal] = secondVarScalarDataInfo_[i].get(volVars);
            }

            // velocity output
            if (velocityOutput.enableOutput())
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    velocityOutput.calculateVelocity(velocity[phaseIdx], elemVolVars, fvGeometry, element, phaseIdx);

            //! the rank
            if (addProcessRank)
                rank[eIdxGlobal] = gridGeom_.gridView().comm().rank();
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
            if (isBox && dim > 1)
            {
                using NestedFunction = VtkNestedFunction<GridView, VertexMapper, std::vector<GlobalPosition>>;
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    sequenceWriter_.addVertexData(std::make_shared<NestedFunction>("velocity_" + std::string(FluidSystem::phaseName(phaseIdx)) + " (m/s)",
                                                                                   gridGeom_.gridView(), gridGeom_.vertexMapper(),
                                                                                   velocity[phaseIdx], dim, dimWorld));
            }
            // cell-centered models
            else
            {
                using NestedFunction = VtkNestedFunction<GridView, ElementMapper, std::vector<GlobalPosition>>;
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    sequenceWriter_.addCellData(std::make_shared<NestedFunction>("velocity_" + std::string(FluidSystem::phaseName(phaseIdx)) + " (m/s)",
                                                                                 gridGeom_.gridView(), gridGeom_.elementMapper(),
                                                                                 velocity[phaseIdx], 0, dimWorld));
            }
        }

        // the process rank
        if (addProcessRank)
            sequenceWriter_.addCellData(rank, "process rank");

        // also register additional (non-standardized) user fields
        for (auto&& field : scalarFields_)
        {
            if (field.first.size() == std::size_t(gridGeom_.gridView().size(0)))
                sequenceWriter_.addCellData(field.first, field.second);
            else if (field.first.size() == std::size_t(gridGeom_.gridView().size(dim)))
                sequenceWriter_.addVertexData(field.first, field.second);
            else
                DUNE_THROW(Dune::RangeError, "Cannot add wrongly sized vtk scalar field!");
        }

        for (auto&& field : vectorFields_)
        {
            if (field.first.size() == std::size_t(gridGeom_.gridView().size(0)))
            {
                using NestedFunction = VtkNestedFunction<GridView, ElementMapper, std::vector<GlobalPosition>>;
                sequenceWriter_.addCellData(std::make_shared<NestedFunction>(field.second,
                                                                             gridGeom_.gridView(), gridGeom_.elementMapper(),
                                                                             field.first, 0, dimWorld));
            }
            else if (field.first.size() == std::size_t(gridGeom_.gridView().size(dim)))
            {
                using NestedFunction = VtkNestedFunction<GridView, VertexMapper, std::vector<GlobalPosition>>;
                sequenceWriter_.addVertexData(std::make_shared<NestedFunction>(field.second,
                                                                               gridGeom_.gridView(), gridGeom_.vertexMapper(),
                                                                               field.first, dim, dimWorld));
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

        //! output
        timer.stop();
        if (verbose_)
        {
            std::cout << "Writing output for problem \"" << name_ << "\". Took " << timer.elapsed() << " seconds." << std::endl;
        }
    }

    //! clear all data in the writer
    void clear()
    {
        writer_->clear();
        scalarFields_.clear();
        vectorFields_.clear();
    }

    /*!
     * \brief This method writes the complete state of the vtk writer
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
     * \brief This method restores the complete state of the vtk writer
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

    //! return the number of dofs
    std::size_t numDofs_() const
    {
        return gridGeom_.gridView().size(dofCodim);
    }

     /*!
     * \brief Helper function to retrieve privar data from the solution vector
     *        May be specialized.
     *
     * \param dofIdxGlobal The global dof index
     * \param pvIdx The index of the primary variable
     */
    auto getPriVarData_(const std::size_t dofIdxGlobal, const std::size_t pvIdx)
    {
        return sol_[dofIdxGlobal][pvIdx];
    }

    //! Returns the implementation of the output module (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    const Problem& problem_;
    const FVGridGeometry& gridGeom_;
    const GridVariables& gridVariables_;
    const SolutionVector& sol_;

    std::string name_;
    bool verbose_;

    std::shared_ptr<Dune::VTKWriter<GridView>> writer_;
    Dune::VTKSequenceWriter<GridView> sequenceWriter_;

    std::vector<PriVarScalarDataInfo> priVarScalarDataInfo_;
    std::vector<PriVarVectorDataInfo> priVarVectorDataInfo_;
    std::vector<SecondVarScalarDataInfo> secondVarScalarDataInfo_;

    std::list<std::pair<std::vector<Scalar>, std::string>> scalarFields_;
    std::list<std::pair<std::vector<GlobalPosition>, std::string>> vectorFields_;
};

} // end namespace Dumux

#endif
