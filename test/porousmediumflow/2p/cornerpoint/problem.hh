// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \ingroup TwoPTests
 * \brief The incompressible 2p cornerpoint test problem.
 */

#ifndef DUMUX_TWOP_CORNERPOINT_TEST_PROBLEM_HH
#define DUMUX_TWOP_CORNERPOINT_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief The incompressible 2p cornerpoint test problem.
 */
template<class TypeTag>
class TwoPCornerPointTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum { dimWorld = GridView::dimensionworld };

public:
    TwoPCornerPointTestProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                               std::shared_ptr<typename ParentType::SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        gravity_ = {0, 0, 9.81};
        injectionElement_ = getParam<int>("Problem.InjectionElement");
        injectionRate_ = getParam<Scalar>("Problem.InjectionRate");
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub-control volume face
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes bcTypes;

        // set no-flux on top and bottom, hydrostatic on the rest
        // use the scvf normal to decide
        const auto& normal = scvf.unitOuterNormal();
        using std::abs;
        if (abs(normal[dimWorld-1]) > 0.5)
            bcTypes.setAllNeumann();
        else
            bcTypes.setAllDirichlet();

        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initialAtPos(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    //! \copydoc Dumux::FVProblem::source()
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        NumEqVector values(0.0);

        int eIdx = this->gridGeometry().gridView().indexSet().index(element);
        if (eIdx == injectionElement_)
            values[FluidSystem::phase1Idx] = injectionRate_/element.geometry().volume();

        return values;
    }

    const GlobalPosition& gravity() const
    {
        return gravity_;
    }

    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        // hydrostatic pressure
        Scalar densityW = 1000;
        values[Indices::pressureIdx] = 1e5 + densityW*(this->spatialParams().gravity(globalPos)*globalPos);
        values[Indices::saturationIdx] = 0.0;

        return values;
    }



    /*!
     * \brief Appends all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    template<class VTKWriter>
    void addFieldsToWriter(VTKWriter& vtk)
    {
        const auto numElements = this->gridGeometry().gridView().size(0);

        permX_.resize(numElements);
        permZ_.resize(numElements);

        vtk.addField(permX_, "PERMX [mD]");
        vtk.addField(permZ_, "PERMZ [mD]");

        const auto& gridView = this->gridGeometry().gridView();
        for (const auto& element : elements(gridView))
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);

            // transfer output to mD = 9.86923e-16 m^2
            permX_[eIdx] = this->spatialParams().permeabilityX(eIdx)/9.86923e-16;
            permZ_[eIdx] = this->spatialParams().permeabilityZ(eIdx)/9.86923e-16;
        }
    }

private:
    GlobalPosition gravity_;
    int injectionElement_;
    Scalar injectionRate_;
    std::vector<Scalar> permX_, permZ_;
};

} // end namespace Dumux


#endif
