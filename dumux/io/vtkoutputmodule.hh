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
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 */
#ifndef DUMUX_IO_VTK_OUTPUT_MODULE_HH
#define DUMUX_IO_VTK_OUTPUT_MODULE_HH

#include <functional>
#include <memory>
#include <string>

#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/common/partitionset.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/io/format.hh>
#include <dumux/discretization/method.hh>

#include "vtkfunction.hh"
#include "velocityoutput.hh"

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 * \note This is a base class providing only rudimentary features
 */
template<class GridGeometry>
class VtkOutputModuleBase
{
    using GridView = typename GridGeometry::GridView;
    using Field = Vtk::template Field<GridView>;
    static constexpr int dim = GridView::dimension;

public:
    //! export field type
    enum class FieldType : unsigned int
    {
        element, vertex, automatic
    };

    VtkOutputModuleBase(const GridGeometry& gridGeometry,
                        const std::string& name,
                        const std::string& paramGroup = "",
                        Dune::VTK::DataMode dm = Dune::VTK::conforming,
                        bool verbose = true)
    : gridGeometry_(gridGeometry)
    , name_(name)
    , paramGroup_(paramGroup)
    , dm_(dm)
    , verbose_(gridGeometry.gridView().comm().rank() == 0 && verbose)
    {
        const auto precisionString = getParamFromGroup<std::string>(paramGroup, "Vtk.Precision", "Float32");
        precision_ = Dumux::Vtk::stringToPrecision(precisionString);
        const auto coordPrecision = Dumux::Vtk::stringToPrecision(getParamFromGroup<std::string>(paramGroup, "Vtk.CoordPrecision", precisionString));
        writer_ = std::make_shared<Dune::VTKWriter<GridView>>(gridGeometry.gridView(), dm, coordPrecision);
        sequenceWriter_ = std::make_unique<Dune::VTKSequenceWriter<GridView>>(writer_, name);
    }

    virtual ~VtkOutputModuleBase() = default;

    //! the parameter group for getting parameter from the parameter tree
    const std::string& paramGroup() const
    { return paramGroup_; }

    /*!
     * \brief Add a scalar or vector valued vtk field
     *
     * \param v The field to be added. Can be any indexable container. Its value type can be a number or itself an indexable container.
     * \param name The name of the field
     * \param fieldType The type of the field.
     *        This determines whether the values are associated with vertices or elements.
     *        By default, the method automatically deduces the correct type for the given input.
     */
    template<typename Vector>
    void addField(const Vector& v,
                  const std::string& name,
                  FieldType fieldType = FieldType::automatic)
    { addField(v, name, this->precision(), fieldType); }

    /*!
     * \brief Add a scalar or vector valued vtk field
     *
     * \param v The field to be added. Can be any indexable container. Its value type can be a number or itself an indexable container.
     * \param name The name of the field
     * \param fieldType The type of the field.
     *        This determines whether the values are associated with vertices or elements.
     *        By default, the method automatically deduces the correct type for the given input.
     * \param precision The output precision of this field (see Dune::VTK::Precision)
     */
    template<typename Vector>
    void addField(const Vector& v,
                  const std::string& name,
                  Dumux::Vtk::Precision precision,
                  FieldType fieldType = FieldType::automatic)
    {
        // Deduce the number of components from the given vector type
        const auto nComp = getNumberOfComponents_(v);

        const auto numElemDofs = gridGeometry().elementMapper().size();
        const auto numVertexDofs = gridGeometry().vertexMapper().size();

        // Automatically deduce the field type ...
        if(fieldType == FieldType::automatic)
        {
            if(numElemDofs == numVertexDofs)
                DUNE_THROW(Dune::InvalidStateException, "Automatic deduction of FieldType failed. Please explicitly specify FieldType::element or FieldType::vertex.");

            if(v.size() == numElemDofs)
                fieldType = FieldType::element;
            else if(v.size() == numVertexDofs)
                fieldType = FieldType::vertex;
            else
                DUNE_THROW(Dune::RangeError, "Size mismatch of added field!");
        }
        // ... or check if the user-specified type matches the size of v
        else
        {
            if(fieldType == FieldType::element)
                if(v.size() != numElemDofs)
                    DUNE_THROW(Dune::RangeError, "Size mismatch of added field!");

            if(fieldType == FieldType::vertex)
                if(v.size() != numVertexDofs)
                    DUNE_THROW(Dune::RangeError, "Size mismatch of added field!");
        }

        // add the appropriate field
        if (fieldType == FieldType::element)
            fields_.emplace_back(gridGeometry_.gridView(), gridGeometry_.elementMapper(), v, name, nComp, 0, dm_, precision);
        else
            fields_.emplace_back(gridGeometry_.gridView(), gridGeometry_.vertexMapper(), v, name, nComp, dim, dm_, precision);
    }

    //! Write the data for this timestep to file in four steps
    //! (1) We assemble all registered variable fields
    //! (2) We register them with the vtk writer
    //! (3) The writer writes the output for us
    //! (4) Clear the writer for the next time step
    void write(double time, Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
        Dune::Timer timer;

        // write to file depending on data mode
        if (dm_ == Dune::VTK::conforming)
            writeConforming_(time, type);
        else if (dm_ == Dune::VTK::nonconforming)
            writeNonConforming_(time, type);
        else
            DUNE_THROW(Dune::NotImplemented, "Output for provided VTK data mode");

        //! output
        timer.stop();
        if (verbose_)
            std::cout << Fmt::format("Writing output for problem \"{}\". Took {:.2g} seconds.\n", name_, timer.elapsed());
    }

protected:
    const GridGeometry& gridGeometry() const { return gridGeometry_; }

    bool verbose() const { return verbose_; }
    const std::string& name() const { return name_; }
    Dune::VTK::DataMode dataMode() const { return dm_; }
    Dumux::Vtk::Precision precision() const { return precision_; }

    Dune::VTKWriter<GridView>& writer() { return *writer_; }
    Dune::VTKSequenceWriter<GridView>& sequenceWriter() { return *sequenceWriter_; }

    const std::vector<Field>& fields() const { return fields_; }

private:
    //! Assembles the fields and adds them to the writer (conforming output)
    virtual void writeConforming_(double time, Dune::VTK::OutputType type)
    {
        //////////////////////////////////////////////////////////////
        //! (1) Assemble all variable fields and add to writer
        //////////////////////////////////////////////////////////////

        // process rank
        static bool addProcessRank = getParamFromGroup<bool>(this->paramGroup(), "Vtk.AddProcessRank");
        std::vector<int> rank;

       //! Abort if no data was registered
        if (!fields_.empty() || addProcessRank)
        {
            const auto numCells = gridGeometry_.gridView().size(0);

            // maybe allocate space for the process rank
            if (addProcessRank)
            {
                rank.resize(numCells);

                for (const auto& element : elements(gridGeometry_.gridView(), Dune::Partitions::interior))
                {
                    const auto eIdxGlobal = gridGeometry_.elementMapper().index(element);
                    rank[eIdxGlobal] = gridGeometry_.gridView().comm().rank();
                }
            }

            //////////////////////////////////////////////////////////////
            //! (2) Register data fields with the vtk writer
            //////////////////////////////////////////////////////////////

            // the process rank
            if (addProcessRank)
                this->sequenceWriter().addCellData(Field(gridGeometry_.gridView(), gridGeometry_.elementMapper(), rank, "process rank", 1, 0).get());

            // also register additional (non-standardized) user fields if any
            for (auto&& field : fields_)
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

        //////////////////////////////////////////////////////////////
        //! (3) Clear the writer
        //////////////////////////////////////////////////////////////
        this->writer().clear();
    }

    //! Assembles the fields and adds them to the writer (nonconforming output)
    virtual void writeNonConforming_(double time, Dune::VTK::OutputType type)
    {
        DUNE_THROW(Dune::NotImplemented, "Non-conforming VTK output");
    }

    //! Deduces the number of components of the value type of a vector of values
    template<class Vector>
    std::size_t getNumberOfComponents_(const Vector& v)
    {
        if constexpr (IsIndexable<decltype(std::declval<Vector>()[0])>::value)
            return v[0].size();
        else
            return 1;
    }

    const GridGeometry& gridGeometry_;
    std::string name_;
    const std::string paramGroup_;
    Dune::VTK::DataMode dm_;
    bool verbose_;
    Dumux::Vtk::Precision precision_;

    std::shared_ptr<Dune::VTKWriter<GridView>> writer_;
    std::unique_ptr<Dune::VTKSequenceWriter<GridView>> sequenceWriter_;

    std::vector<Field> fields_; //!< Registered scalar and vector fields
};

/*!
 * \ingroup InputOutput
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 *
 * \tparam GridVariables The grid variables
 * \tparam SolutionVector The solution vector
 *
 * Handles the output of scalar and vector fields to VTK formatted file for multiple
 * variables and timesteps. Certain predefined fields can be registered on
 * initialization and/or be turned on/off using the designated properties. Additionally
 * non-standardized scalar and vector fields can be added to the writer manually.
 */
template<class GridVariables, class SolutionVector>
class VtkOutputModule : public VtkOutputModuleBase<typename GridVariables::GridGeometry>
{
    using ParentType = VtkOutputModuleBase<typename GridVariables::GridGeometry>;
    using GridGeometry = typename GridVariables::GridGeometry;

    using VV = typename GridVariables::VolumeVariables;
    using Scalar = typename GridVariables::Scalar;

    using GridView = typename GridGeometry::GridView;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using VolVarsVector = Dune::FieldVector<Scalar, dimWorld>;

    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethod::box;
    static constexpr int dofCodim = isBox ? dim : 0;

    struct VolVarScalarDataInfo { std::function<Scalar(const VV&)> get; std::string name; Dumux::Vtk::Precision precision_; };
    struct VolVarVectorDataInfo { std::function<VolVarsVector(const VV&)> get; std::string name; Dumux::Vtk::Precision precision_; };
    using Field = Vtk::template Field<GridView>;

    using VelocityOutputType = Dumux::VelocityOutput<GridVariables>;

public:
    //! export type of the volume variables for the outputfields
    using VolumeVariables = VV;

    VtkOutputModule(const GridVariables& gridVariables,
                    const SolutionVector& sol,
                    const std::string& name,
                    const std::string& paramGroup = "",
                    Dune::VTK::DataMode dm = Dune::VTK::conforming,
                    bool verbose = true)
    : ParentType(gridVariables.gridGeometry(), name, paramGroup, dm, verbose)
    , gridVariables_(gridVariables)
    , sol_(sol)
    , velocityOutput_(std::make_shared<VelocityOutputType>())
    {}

    //////////////////////////////////////////////////////////////////////////////////////////////
    //! Methods to conveniently add primary and secondary variables upon initialization
    //! Do not call these methods after initialization i.e. _not_ within the time loop
    //////////////////////////////////////////////////////////////////////////////////////////////

    /*!
     * \brief Add a velocity output policy
     *
     * \param velocityOutput the output policy
     * \note the default policy does not add any velocity output
     */
    void addVelocityOutput(std::shared_ptr<VelocityOutputType> velocityOutput)
    { velocityOutput_ = velocityOutput; }

    //! Output a scalar volume variable
    //! \param name The name of the vtk field
    //! \param f A function taking a VolumeVariables object and returning the desired scalar
    void addVolumeVariable(std::function<Scalar(const VolumeVariables&)>&& f,
                                                const std::string& name)
    {
        volVarScalarDataInfo_.push_back(VolVarScalarDataInfo{f, name, this->precision()});
    }

    //! Add a vector-valued variable
    //! \param f A function taking a VolumeVariables object and returning the desired vector
    //! \param name The name of the vtk field
    //! \note This method is only available for dimWorld > 1. For 1-D problems, the overload for volVar methods returning a Scalar will be used.
    template<class VVV = VolVarsVector, typename std::enable_if_t<(VVV::dimension > 1), int> = 0>
    void addVolumeVariable(std::function<VolVarsVector(const VolumeVariables&)>&& f,
                                                       const std::string& name)
    {
        volVarVectorDataInfo_.push_back(VolVarVectorDataInfo{f, name, this->precision()});
    }

protected:
    // some return functions for differing implementations to use
    const auto& problem() const { return gridVariables_.curGridVolVars().problem(); }
    const GridVariables& gridVariables() const { return gridVariables_; }
    const GridGeometry& gridGeometry() const { return gridVariables_.gridGeometry(); }
    const SolutionVector& sol() const { return sol_; }

    const std::vector<VolVarScalarDataInfo>& volVarScalarDataInfo() const { return volVarScalarDataInfo_; }
    const std::vector<VolVarVectorDataInfo>& volVarVectorDataInfo() const { return volVarVectorDataInfo_; }

    using VelocityOutput = VelocityOutputType;
    const VelocityOutput& velocityOutput() const { return *velocityOutput_; }

private:

    //! Assembles the fields and adds them to the writer (conforming output)
    void writeConforming_(double time, Dune::VTK::OutputType type) override
    {
        const Dune::VTK::DataMode dm = Dune::VTK::conforming;
        //////////////////////////////////////////////////////////////
        //! (1) Assemble all variable fields and add to writer
        //////////////////////////////////////////////////////////////

        // instantiate the velocity output
        using VelocityVector = typename VelocityOutput::VelocityVector;
        std::vector<VelocityVector> velocity(velocityOutput_->numFluidPhases());

        // process rank
        static bool addProcessRank = getParamFromGroup<bool>(this->paramGroup(), "Vtk.AddProcessRank");
        std::vector<double> rank;

        // volume variable data
        std::vector<std::vector<Scalar>> volVarScalarData;
        std::vector<std::vector<VolVarsVector>> volVarVectorData;

        //! Abort if no data was registered
        if (!volVarScalarDataInfo_.empty()
            || !volVarVectorDataInfo_.empty()
            || !this->fields().empty()
            || velocityOutput_->enableOutput()
            || addProcessRank)
        {
            const auto numCells = gridGeometry().gridView().size(0);
            const auto numDofs = numDofs_();

            // get fields for all volume variables
            if (!volVarScalarDataInfo_.empty())
                volVarScalarData.resize(volVarScalarDataInfo_.size(), std::vector<Scalar>(numDofs));
            if (!volVarVectorDataInfo_.empty())
                volVarVectorData.resize(volVarVectorDataInfo_.size(), std::vector<VolVarsVector>(numDofs));

            if (velocityOutput_->enableOutput())
            {
                for (int phaseIdx = 0; phaseIdx < velocityOutput_->numFluidPhases(); ++phaseIdx)
                {
                    if(isBox && dim == 1)
                        velocity[phaseIdx].resize(numCells);
                    else
                        velocity[phaseIdx].resize(numDofs);
                }
            }

            // maybe allocate space for the process rank
            if (addProcessRank) rank.resize(numCells);

            for (const auto& element : elements(gridGeometry().gridView(), Dune::Partitions::interior))
            {
                const auto eIdxGlobal = gridGeometry().elementMapper().index(element);

                auto fvGeometry = localView(gridGeometry());
                auto elemVolVars = localView(gridVariables_.curGridVolVars());

                // If velocity output is enabled we need to bind to the whole stencil
                // otherwise element-local data is sufficient
                if (velocityOutput_->enableOutput())
                {
                    fvGeometry.bind(element);
                    elemVolVars.bind(element, fvGeometry, sol_);
                }
                else
                {
                    fvGeometry.bindElement(element);
                    elemVolVars.bindElement(element, fvGeometry, sol_);
                }

                if (!volVarScalarDataInfo_.empty()
                    || !volVarVectorDataInfo_.empty())
                {
                    for (auto&& scv : scvs(fvGeometry))
                    {
                        const auto dofIdxGlobal = scv.dofIndex();
                        const auto& volVars = elemVolVars[scv];

                        // get the scalar-valued data
                        for (std::size_t i = 0; i < volVarScalarDataInfo_.size(); ++i)
                            volVarScalarData[i][dofIdxGlobal] = volVarScalarDataInfo_[i].get(volVars);

                        // get the vector-valued data
                        for (std::size_t i = 0; i < volVarVectorDataInfo_.size(); ++i)
                            volVarVectorData[i][dofIdxGlobal] = volVarVectorDataInfo_[i].get(volVars);
                    }
                }

                // velocity output
                if (velocityOutput_->enableOutput())
                {
                    auto elemFluxVarsCache = localView(gridVariables_.gridFluxVarsCache());
                    elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

                    for (int phaseIdx = 0; phaseIdx < velocityOutput_->numFluidPhases(); ++phaseIdx)
                        velocityOutput_->calculateVelocity(velocity[phaseIdx], element, fvGeometry, elemVolVars, elemFluxVarsCache, phaseIdx);
                }

                //! the rank
                if (addProcessRank)
                    rank[eIdxGlobal] = static_cast<double>(gridGeometry().gridView().comm().rank());
            }

            //////////////////////////////////////////////////////////////
            //! (2) Register data fields with the vtk writer
            //////////////////////////////////////////////////////////////

            // volume variables if any
            if (isBox)
            {
                for (std::size_t i = 0; i < volVarScalarDataInfo_.size(); ++i)
                    this->sequenceWriter().addVertexData( Field(gridGeometry().gridView(), gridGeometry().vertexMapper(), volVarScalarData[i],
                                                         volVarScalarDataInfo_[i].name, /*numComp*/1, /*codim*/dim, dm, this->precision()).get() );
                for (std::size_t i = 0; i < volVarVectorDataInfo_.size(); ++i)
                    this->sequenceWriter().addVertexData( Field(gridGeometry().gridView(), gridGeometry().vertexMapper(), volVarVectorData[i],
                                                         volVarVectorDataInfo_[i].name, /*numComp*/dimWorld, /*codim*/dim, dm, this->precision()).get() );
            }
            else
            {
                for (std::size_t i = 0; i < volVarScalarDataInfo_.size(); ++i)
                    this->sequenceWriter().addCellData( Field(gridGeometry().gridView(), gridGeometry().elementMapper(), volVarScalarData[i],
                                                       volVarScalarDataInfo_[i].name, /*numComp*/1, /*codim*/0,dm, this->precision()).get() );
                for (std::size_t i = 0; i < volVarVectorDataInfo_.size(); ++i)
                    this->sequenceWriter().addCellData( Field(gridGeometry().gridView(), gridGeometry().elementMapper(), volVarVectorData[i],
                                                       volVarVectorDataInfo_[i].name, /*numComp*/dimWorld, /*codim*/0,dm, this->precision()).get() );
            }

            // the velocity field
            if (velocityOutput_->enableOutput())
            {
                if (isBox && dim > 1)
                {
                    for (int phaseIdx = 0; phaseIdx < velocityOutput_->numFluidPhases(); ++phaseIdx)
                        this->sequenceWriter().addVertexData( Field(gridGeometry().gridView(), gridGeometry().vertexMapper(), velocity[phaseIdx],
                                                             "velocity_" + velocityOutput_->phaseName(phaseIdx) + " (m/s)",
                                                             /*numComp*/dimWorld, /*codim*/dim, dm, this->precision()).get() );
                }
                // cell-centered models
                else
                {
                    for (int phaseIdx = 0; phaseIdx < velocityOutput_->numFluidPhases(); ++phaseIdx)
                        this->sequenceWriter().addCellData( Field(gridGeometry().gridView(), gridGeometry().elementMapper(), velocity[phaseIdx],
                                                           "velocity_" + velocityOutput_->phaseName(phaseIdx) + " (m/s)",
                                                           /*numComp*/dimWorld, /*codim*/0, dm, this->precision()).get() );
                }
            }

            // the process rank
            if (addProcessRank)
                this->sequenceWriter().addCellData(Field(gridGeometry().gridView(), gridGeometry().elementMapper(), rank, "process rank", 1, 0).get());

            // also register additional (non-standardized) user fields if any
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

        //////////////////////////////////////////////////////////////
        //! (3) Clear the writer
        //////////////////////////////////////////////////////////////
        this->writer().clear();
    }

    //! Assembles the fields and adds them to the writer (nonconforming output)
    void writeNonConforming_(double time, Dune::VTK::OutputType type) override
    {
        const Dune::VTK::DataMode dm = Dune::VTK::nonconforming;

        if(!isBox)
            DUNE_THROW(Dune::InvalidStateException, "Non-conforming output makes no sense for cell-centered schemes!");

        //////////////////////////////////////////////////////////////
        //! (1) Assemble all variable fields and add to writer
        //////////////////////////////////////////////////////////////

        // check the velocity output
        bool enableVelocityOutput = getParamFromGroup<bool>(this->paramGroup(), "Vtk.AddVelocity");
        if (enableVelocityOutput == true && !velocityOutput_->enableOutput())
            std::cerr << "Warning! Velocity output was enabled in the input file"
                      << " but no velocity output policy was set for the VTK output module:"
                      << " There will be no velocity output."
                      << " Use the addVelocityOutput member function of the VTK output module." << std::endl;
        using VelocityVector = typename VelocityOutput::VelocityVector;
        std::vector<VelocityVector> velocity(velocityOutput_->numFluidPhases());

        // process rank
        static bool addProcessRank = getParamFromGroup<bool>(this->paramGroup(), "Vtk.AddProcessRank");
        std::vector<double> rank;

        // volume variable data (indexing: volvardata/element/localcorner)
        using ScalarDataContainer = std::vector< std::vector<Scalar> >;
        using VectorDataContainer = std::vector< std::vector<VolVarsVector> >;
        std::vector< ScalarDataContainer > volVarScalarData;
        std::vector< VectorDataContainer > volVarVectorData;

        //! Abort if no data was registered
        if (!volVarScalarDataInfo_.empty()
            || !volVarVectorDataInfo_.empty()
            || !this->fields().empty()
            || velocityOutput_->enableOutput()
            || addProcessRank)
        {
            const auto numCells = gridGeometry().gridView().size(0);
            const auto numDofs = numDofs_();

            // get fields for all volume variables
            if (!volVarScalarDataInfo_.empty())
                volVarScalarData.resize(volVarScalarDataInfo_.size(), ScalarDataContainer(numCells));
            if (!volVarVectorDataInfo_.empty())
                volVarVectorData.resize(volVarVectorDataInfo_.size(), VectorDataContainer(numCells));

            if (velocityOutput_->enableOutput())
            {
                for (int phaseIdx = 0; phaseIdx < velocityOutput_->numFluidPhases(); ++phaseIdx)
                {
                    if(isBox && dim == 1)
                        velocity[phaseIdx].resize(numCells);
                    else
                        velocity[phaseIdx].resize(numDofs);
                }
            }

            // maybe allocate space for the process rank
            if (addProcessRank) rank.resize(numCells);

            for (const auto& element : elements(gridGeometry().gridView(), Dune::Partitions::interior))
            {
                const auto eIdxGlobal = gridGeometry().elementMapper().index(element);
                const auto numCorners = element.subEntities(dim);

                auto fvGeometry = localView(gridGeometry());
                auto elemVolVars = localView(gridVariables_.curGridVolVars());

                // resize element-local data containers
                for (std::size_t i = 0; i < volVarScalarDataInfo_.size(); ++i)
                    volVarScalarData[i][eIdxGlobal].resize(numCorners);
                for (std::size_t i = 0; i < volVarVectorDataInfo_.size(); ++i)
                    volVarVectorData[i][eIdxGlobal].resize(numCorners);

                // If velocity output is enabled we need to bind to the whole stencil
                // otherwise element-local data is sufficient
                if (velocityOutput_->enableOutput())
                {
                    fvGeometry.bind(element);
                    elemVolVars.bind(element, fvGeometry, sol_);
                }
                else
                {
                    fvGeometry.bindElement(element);
                    elemVolVars.bindElement(element, fvGeometry, sol_);
                }

                if (!volVarScalarDataInfo_.empty()
                    || !volVarVectorDataInfo_.empty())
                {
                    for (auto&& scv : scvs(fvGeometry))
                    {
                        const auto& volVars = elemVolVars[scv];

                        // get the scalar-valued data
                        for (std::size_t i = 0; i < volVarScalarDataInfo_.size(); ++i)
                            volVarScalarData[i][scv.elementIndex()][scv.localDofIndex()] = volVarScalarDataInfo_[i].get(volVars);

                        // get the vector-valued data
                        for (std::size_t i = 0; i < volVarVectorDataInfo_.size(); ++i)
                            volVarVectorData[i][scv.elementIndex()][scv.localDofIndex()] = volVarVectorDataInfo_[i].get(volVars);
                    }
                }

                // velocity output
                if (velocityOutput_->enableOutput())
                {
                    auto elemFluxVarsCache = localView(gridVariables_.gridFluxVarsCache());
                    elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

                    for (int phaseIdx = 0; phaseIdx < velocityOutput_->numFluidPhases(); ++phaseIdx)
                        velocityOutput_->calculateVelocity(velocity[phaseIdx], element, fvGeometry, elemVolVars, elemFluxVarsCache, phaseIdx);
                }

                //! the rank
                if (addProcessRank)
                    rank[eIdxGlobal] = static_cast<double>(gridGeometry().gridView().comm().rank());
            }

            //////////////////////////////////////////////////////////////
            //! Register data fields with the vtk writer
            //////////////////////////////////////////////////////////////

            // volume variables if any
            for (std::size_t i = 0; i < volVarScalarDataInfo_.size(); ++i)
                this->sequenceWriter().addVertexData( Field(gridGeometry().gridView(), gridGeometry().elementMapper(), volVarScalarData[i],
                                                     volVarScalarDataInfo_[i].name, /*numComp*/1, /*codim*/dim, /*nonconforming*/dm, this->precision()).get() );

            for (std::size_t i = 0; i < volVarVectorDataInfo_.size(); ++i)
                this->sequenceWriter().addVertexData( Field(gridGeometry().gridView(), gridGeometry().elementMapper(), volVarVectorData[i],
                                                     volVarVectorDataInfo_[i].name, /*numComp*/dimWorld, /*codim*/dim, /*nonconforming*/dm, this->precision()).get() );

            // the velocity field
            if (velocityOutput_->enableOutput())
            {
                // node-wise velocities
                if (dim > 1)
                    for (int phaseIdx = 0; phaseIdx < velocityOutput_->numFluidPhases(); ++phaseIdx)
                        this->sequenceWriter().addVertexData( Field(gridGeometry().gridView(), gridGeometry().vertexMapper(), velocity[phaseIdx],
                                                             "velocity_" + velocityOutput_->phaseName(phaseIdx) + " (m/s)",
                                                             /*numComp*/dimWorld, /*codim*/dim, dm, this->precision()).get() );

                // cell-wise velocities
                else
                    for (int phaseIdx = 0; phaseIdx < velocityOutput_->numFluidPhases(); ++phaseIdx)
                        this->sequenceWriter().addCellData( Field(gridGeometry().gridView(), gridGeometry().elementMapper(), velocity[phaseIdx],
                                                           "velocity_" + velocityOutput_->phaseName(phaseIdx) + " (m/s)",
                                                           /*numComp*/dimWorld, /*codim*/0,dm, this->precision()).get());
            }

            // the process rank
            if (addProcessRank)
                this->sequenceWriter().addCellData( Field(gridGeometry().gridView(), gridGeometry().elementMapper(), rank, "process rank", 1, 0).get() );

            // also register additional (non-standardized) user fields if any
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

        //////////////////////////////////////////////////////////////
        //! (3) Clear the writer
        //////////////////////////////////////////////////////////////
        this->writer().clear();
    }

    //! return the number of dofs, we only support vertex and cell data
    std::size_t numDofs_() const { return dofCodim == dim ? gridGeometry().vertexMapper().size() : gridGeometry().elementMapper().size(); }

    const GridVariables& gridVariables_;
    const SolutionVector& sol_;

    std::vector<VolVarScalarDataInfo> volVarScalarDataInfo_; //!< Registered volume variables (scalar)
    std::vector<VolVarVectorDataInfo> volVarVectorDataInfo_; //!< Registered volume variables (vector)

    std::shared_ptr<VelocityOutput> velocityOutput_; //!< The velocity output policy
};

} // end namespace Dumux

#endif
