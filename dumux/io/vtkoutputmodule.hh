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
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 */
#ifndef VTK_OUTPUT_MODULE_HH
#define VTK_OUTPUT_MODULE_HH

#include <functional>

#include <dune/common/version.hh>
#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/io/file/vtk/function.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/common/partitionset.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {
namespace Vtk {

//! struct that can hold any field that fulfills the VTKFunction interface
template<class GridView>
class Field
{
    enum { dim = GridView::dimension };
    using ctype = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;

    // a VTK function that supports both scalar and vector values for each element
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
    template <typename F>
    struct VectorP0VTKFunction : Dune::VTKFunction<GridView>
    {
        using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
#else
    template<typename F, typename Mapper>
    struct VectorP0VTKFunction : Dune::VTKFunction<GridView>
    {
#endif
    public:
        //! return number of components
        virtual int ncomps () const
        { return nComps_; }

        //! evaluate
        virtual double evaluate (int mycomp, const Element& e,
                                 const Dune::FieldVector<ctype,dim>&) const
        { return accessChooser_(mycomp, mapper_.index(e), Dune::is_indexable<decltype(field_[0])>()); }

        //! get name
        virtual std::string name () const
        { return name_; }

        VectorP0VTKFunction(const GridView &gridView, const Mapper& mapper, const F& field, const std::string& name, int nComps)
        : field_(field), name_(name), nComps_(nComps), mapper_(mapper)
        {
            if (field.size()!=(unsigned int)(gridView.size(0)))
                DUNE_THROW(Dune::IOError, "NestedP0VTKFunction: size mismatch");
        }
    private:
        double accessChooser_(int mycomp, int i, std::true_type) const
        { return field_[i][mycomp]; }

        double accessChooser_(int mycomp, int i, std::false_type) const
        { return field_[i]; }

        const F& field_;
        const std::string name_;
        int nComps_;
        const Mapper& mapper_;
    };

    // a VTK function that supports both scalar and vector values for each vertex
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
    template <typename F>
    struct VectorP1VTKFunction : Dune::VTKFunction<GridView>
    {
        using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
#else
    template<typename F, typename Mapper>
    struct VectorP1VTKFunction : Dune::VTKFunction<GridView>
    {
#endif
    public:
        //! return number of components
        virtual int ncomps () const
        { return nComps_; }

        //! evaluate
        virtual double evaluate (int mycomp, const Element& e,
                                 const Dune::FieldVector<ctype,dim>& xi) const
        {
            const unsigned int dim = Element::mydimension;
            const unsigned int nVertices = e.subEntities(dim);

            std::vector<Dune::FieldVector<ctype, 1>> cornerValues(nVertices);
            for (unsigned i = 0; i < nVertices; ++i)
                cornerValues[i] = accessChooser_(mycomp, mapper_.subIndex(e, i, dim), Dune::is_indexable<decltype(field_[0])>());

            // (Ab)use the MultiLinearGeometry class to do multi-linear interpolation between scalars
            const Dune::MultiLinearGeometry<ctype, dim, 1> interpolation(e.type(), std::move(cornerValues));
            return interpolation.global(xi);
        }

        //! get name
        virtual std::string name () const
        { return name_; }

        VectorP1VTKFunction(const GridView &gridView, const Mapper& mapper, const F& field, const std::string& name, int nComps)
        : field_(field), name_(name), nComps_(nComps), mapper_(mapper)
        {
            if (field.size()!=(unsigned int)(gridView.size(GridView::dimension)))
                DUNE_THROW(Dune::IOError, "NestedP1VTKFunction: size mismatch");
        }
    private:
        double accessChooser_(int mycomp, int i, std::true_type) const
        { return field_[i][mycomp]; }

        double accessChooser_(int mycomp, int i, std::false_type) const
        { return field_[i]; }

        const F& field_;
        const std::string name_;
        int nComps_;
        const Mapper& mapper_;
    };

public:
    // template constructor selects the right VTKFunction implementation
    template <typename F, class Mapper>
    Field(const GridView& gridView, const Mapper& mapper, F const& f,
          const std::string& name, int numComp = 1, int codim = 0)
    : codim_(codim)
    {
        if (codim == GridView::dimension)
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
            field_ = std::make_shared<VectorP1VTKFunction<F>>(gridView, mapper, f, name, numComp);
#else
            field_ = std::make_shared<VectorP1VTKFunction<F, Mapper>>(gridView, mapper, f, name, numComp);
#endif
        else if (codim == 0)
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
            field_ = std::make_shared<VectorP0VTKFunction<F>>(gridView, mapper, f, name, numComp);
#else
            field_ = std::make_shared<VectorP0VTKFunction<F, Mapper>>(gridView, mapper, f, name, numComp);
#endif
        else
            DUNE_THROW(Dune::NotImplemented, "Only element or vertex quantities allowed.");
    }

    virtual std::string name () const
    { return field_->name(); }

    virtual int ncomps() const
    { return field_->ncomps(); }

    virtual double evaluate(int mycomp,
                            const Element &element,
                            const Dune::FieldVector< ctype, dim > &xi) const
    { return field_->evaluate(mycomp, element, xi); }

    int codim() const
    { return codim_; }

    const std::shared_ptr<Dune::VTKFunction<GridView>>& get() const
    { return field_; }

private:
    int codim_;
    // can point to anything fulfilling the VTKFunction interface
    std::shared_ptr<Dune::VTKFunction<GridView>> field_;
};

} // end namespace Vtk

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
//    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
//    using VelocityOutput = typename GET_PROP_TYPE(TypeTag, VelocityOutput);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using ElementMapper = typename FVGridGeometry::ElementMapper;
    using VertexMapper = typename FVGridGeometry::VertexMapper;
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, DiscretizationMethod) == DiscretizationMethods::Box;
    static constexpr int dofCodim = isBox ? dim : 0;

    struct VolVarScalarDataInfo { std::function<Scalar(const VolumeVariables&)> get; std::string name; };
    using Field = Vtk::template Field<GridView>;

public:

    enum class FieldType : unsigned int
    {
        element, vertex, automatic
    };

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
    //! Do not call these methods after initialization i.e. _not_ within the time loop
    //////////////////////////////////////////////////////////////////////////////////////////////

    //! Output a scalar volume variable
    //! \param name The name of the vtk field
    //! \param f A function taking a VolumeVariables object and returning the desired scalar
    DUNE_DEPRECATED_MSG("Will be removed after next release. Please use addVolumeVariable(function, name) instead!")
    void addSecondaryVariable(const std::string& name, std::function<Scalar(const VolumeVariables&)>&& f)
    {
        volVarScalarDataInfo_.push_back(VolVarScalarDataInfo{f, name});
    }

    //! Output a scalar volume variable
    //! \param name The name of the vtk field
    //! \param f A function taking a VolumeVariables object and returning the desired scalar
    void addVolumeVariable(std::function<Scalar(const VolumeVariables&)>&& f, const std::string& name)
    {
        volVarScalarDataInfo_.push_back(VolVarScalarDataInfo{f, name});
    }

    //! Add a scalar or vector valued vtk field
    //! \param v The field to be added. Can be any indexable container. Its value type can be a number or itself an indexable container.
    //! \param name The name of the field
    //! \param fieldType The type of the field.
    //!        This determines whether the values are associated with vertices or elements.
    //!        By default, the method automatically deduces the correct type for the given input.
    template<typename Vector>
    void addField(const Vector& v, const std::string& name, FieldType fieldType = FieldType::automatic)
    {
        // Deduce the number of components from the given vector type
        const auto nComp = getNumberOfComponents_(v);

        const auto numElements = gridGeom_.gridView().size(0);
        const auto numVertices = gridGeom_.gridView().size(dim);

        // Automatically deduce the field type ...
        if(fieldType == FieldType::automatic)
        {
            if(numElements == numVertices)
                DUNE_THROW(Dune::InvalidStateException, "Automatic deduction of FieldType failed. Please explicitly specify FieldType::element or FieldType::vertex.");

            if(v.size() == numElements)
                fieldType = FieldType::element;
            else if(v.size() == numVertices)
                fieldType = FieldType::vertex;
            else
                DUNE_THROW(Dune::RangeError, "Size mismatch of added field!");
        }
        // ... or check if the user-specified type matches the size of v
        else
        {
            if(fieldType == FieldType::element)
                if(v.size() != numElements)
                    DUNE_THROW(Dune::RangeError, "Size mismatch of added field!");

            if(fieldType == FieldType::vertex)
                if(v.size() != numVertices)
                    DUNE_THROW(Dune::RangeError, "Size mismatch of added field!");
        }

        // add the appropriate field
        if (fieldType == FieldType::element)
            fields_.emplace_back(gridGeom_.gridView(), gridGeom_.elementMapper(), v, name, nComp, 0);
        else
            fields_.emplace_back(gridGeom_.gridView(), gridGeom_.vertexMapper(), v, name, nComp, dim);
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

        //////////////////////////////////////////////////////////////
        //! (1) Assemble all variable fields with registered info
        //////////////////////////////////////////////////////////////

        // instatiate the velocity output
        //VelocityOutput velocityOutput(problem_, gridGeom_, gridVariables_, sol_);
        std::array<std::vector<GlobalPosition>, numPhases> velocity;

        // process rank
        static const std::string modelParamGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        static bool addProcessRank = getParamFromGroup<bool>(modelParamGroup, "Vtk.AddProcessRank");
        std::vector<double> rank;

        // volume variable data
        std::vector<std::vector<Scalar>> volVarScalarData;

        //! Abort if no data was registered
        if (!volVarScalarDataInfo_.empty()
            || !fields_.empty()
//            || velocityOutput.enableOutput()
            || addProcessRank)
        {
            const auto numCells = gridGeom_.gridView().size(0);
            const auto numDofs = numDofs_();

            // get fields for all volume variables
            if (!volVarScalarDataInfo_.empty())
                volVarScalarData.resize(volVarScalarDataInfo_.size(), std::vector<Scalar>(numDofs));

            if (0)//velocityOutput.enableOutput())
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
            if (addProcessRank) rank.resize(numCells);

            for (const auto& element : elements(gridGeom_.gridView(), Dune::Partitions::interior))
            {
                const auto eIdxGlobal = gridGeom_.elementMapper().index(element);

                auto fvGeometry = localView(gridGeom_);
                auto elemVolVars = localView(gridVariables_.curGridVolVars());

                // If velocity output is enabled we need to bind to the whole stencil
                // otherwise element-local data is sufficient
                if (0)//velocityOutput.enableOutput())
                    fvGeometry.bind(element);
                else
                    fvGeometry.bindElement(element);

                // If velocity output is enabled we need to bind to the whole stencil
                // otherwise element-local data is sufficient
                if (0)//velocityOutput.enableOutput())
                    elemVolVars.bind(element, fvGeometry, sol_);
                else if (!volVarScalarDataInfo_.empty())
                    elemVolVars.bindElement(element, fvGeometry, sol_);

                if (!volVarScalarDataInfo_.empty())
                {
                    for (auto&& scv : scvs(fvGeometry))
                    {
                        const auto dofIdxGlobal = scv.dofIndex();
                        const auto& volVars = elemVolVars[scv];

                        for (std::size_t i = 0; i < volVarScalarDataInfo_.size(); ++i)
                            volVarScalarData[i][dofIdxGlobal] = volVarScalarDataInfo_[i].get(volVars);
                    }
                }

                // velocity output
//                 if (velocityOutput.enableOutput())
//                     for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
//                         velocityOutput.calculateVelocity(velocity[phaseIdx], elemVolVars, fvGeometry, element, phaseIdx);

                //! the rank
                if (addProcessRank)
                    rank[eIdxGlobal] = static_cast<double>(gridGeom_.gridView().comm().rank());
            }

            //////////////////////////////////////////////////////////////
            //! (2) Register data fields with the vtk writer
            //////////////////////////////////////////////////////////////

            // volume variables if any
            for (std::size_t i = 0; i < volVarScalarDataInfo_.size(); ++i)
                addDofDataForWriter_(sequenceWriter_, volVarScalarData[i], volVarScalarDataInfo_[i].name);

            // the velocity field
            if (0)//velocityOutput.enableOutput())
            {
                if (isBox && dim > 1)
                {
                    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                        sequenceWriter_.addVertexData(Field(gridGeom_.gridView(), gridGeom_.vertexMapper(), velocity[phaseIdx],
                                                            "velocity_" + std::to_string(phaseIdx) + " (m/s)",
                                                            dimWorld, dim).get());
                }
                // cell-centered models
                else
                {
                    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                        sequenceWriter_.addCellData(Field(gridGeom_.gridView(), gridGeom_.elementMapper(), velocity[phaseIdx],
                                                            "velocity_" + std::to_string(phaseIdx) + " (m/s)",
                                                            dimWorld, 0).get());
                }
            }

            // the process rank
            if (addProcessRank)
                sequenceWriter_.addCellData(Field(gridGeom_.gridView(), gridGeom_.elementMapper(), rank, "process rank", 1, 0).get());

            // also register additional (non-standardized) user fields if any
            for (auto&& field : fields_)
            {
                if (field.codim() == 0)
                    sequenceWriter_.addCellData(field.get());
                else if (field.codim() == dim)
                    sequenceWriter_.addVertexData(field.get());
                else
                    DUNE_THROW(Dune::RangeError, "Cannot add wrongly sized vtk scalar field!");
            }
        }

        //////////////////////////////////////////////////////////////
        //! (3) The writer writes the output for us
        //////////////////////////////////////////////////////////////
        sequenceWriter_.write(time, type);

        //////////////////////////////////////////////////////////////
        //! (4) Clear the writer
        //////////////////////////////////////////////////////////////
        writer_->clear();

        //! output
        timer.stop();
        if (verbose_)
        {
            std::cout << "Writing output for problem \"" << name_ << "\". Took " << timer.elapsed() << " seconds." << std::endl;
        }
    }

private:
    //! Deduces the number of components of the value type of a vector of values
    template<class Vector, typename std::enable_if_t<Dune::is_indexable<decltype(std::declval<Vector>()[0])>::value, int> = 0>
    std::size_t getNumberOfComponents_(const Vector& v)
    { return v[0].size(); }

    //! Deduces the number of components of the value type of a vector of values
    template<class Vector, typename std::enable_if_t<!Dune::is_indexable<decltype(std::declval<Vector>()[0])>::value, int> = 0>
    std::size_t getNumberOfComponents_(const Vector& v)
    { return 1; }

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

    const Problem& problem_;
    const FVGridGeometry& gridGeom_;
    const GridVariables& gridVariables_;
    const SolutionVector& sol_;

    std::string name_;
    bool verbose_;

    std::shared_ptr<Dune::VTKWriter<GridView>> writer_;
    Dune::VTKSequenceWriter<GridView> sequenceWriter_;

    std::vector<VolVarScalarDataInfo> volVarScalarDataInfo_; //!< Registered volume variables
    std::vector<Field> fields_; //!< Registered scalar and vector fields
};

} // end namespace Dumux

#endif
