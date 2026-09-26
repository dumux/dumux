// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVETests
 * \brief Radial injection into a homogeneous confined aquifer, solved with the two-phase VE model.
 */

#ifndef DUMUX_TEST_TWOPVE_RADIAL_INJECTION_PROBLEM_HH
#define DUMUX_TEST_TWOPVE_RADIAL_INJECTION_PROBLEM_HH

#include <cmath>
#include <memory>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2pve/finelevelview.hh>

namespace Dumux {

/*!
 * \ingroup TwoPVETests
 * \brief Radial injection into a homogeneous confined aquifer, solved with the two-phase VE model
 *
 * The grid spans the radial and the vertical coordinate and is rotated about the vertical axis.
 * The injected fluid enters through the well at the inner radius, and the pressure at the bottom
 * of the aquifer is fixed at the outer radius.
 */
template<class TypeTag>
class TwoPVERadialInjectionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    static constexpr int pressureIdx = Indices::pressureIdx;
    static constexpr int saturationIdx = Indices::saturationIdx;
    static constexpr int injectedPhaseEqIdx = Indices::conti0EqIdx + FluidSystem::phase1Idx;

public:
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using FineLevelView = TwoPVEFineLevelView<GridGeometry, Scalar, FluidSystem, Indices, SolutionVector, typename SpatialParams::SpatialParamsFine>;

    TwoPVERadialInjectionProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                 std::shared_ptr<SpatialParams> spatialParams,
                                 std::shared_ptr<const FineLevelView> fineLevelView)
    : ParentType(gridGeometry, spatialParams)
    , fineLevelView_(fineLevelView)
    , injectionRate_(getParam<Scalar>("Problem.InjectionRate"))
    , initialPressure_(getParam<Scalar>("Problem.InitialPressure"))
    {}

    //! Returns the fine-level view of the VE model
    const FineLevelView& fineLevelView() const
    { return *fineLevelView_; }

    //! The volumetric injection rate of the injected fluid \f$\mathrm{[m^3/s]}\f$
    Scalar injectionRate() const
    { return injectionRate_; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (onOuterBoundary_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }

    /*!
     * \brief The mass flux of the injected fluid through the well, distributed over the aquifer height
     *
     * \param globalPos the center of the boundary face
     */
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);
        if (onWell_(globalPos))
        {
            const Scalar wellRadius = this->gridGeometry().bBoxMin()[0];
            const Scalar height = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];
            values[injectedPhaseEqIdx] = -injectedDensity_()*injectionRate_/(2.0*M_PI*wellRadius*height);
        }
        return values;
    }

    //! The pressure of the resident fluid at the bottom of the aquifer and a fully water-saturated aquifer
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] = initialPressure_;
        values[saturationIdx] = 0.0;
        return values;
    }

private:
    static constexpr Scalar eps_ = 1e-6;

    bool onWell_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onOuterBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    Scalar injectedDensity_() const
    {
        typename FluidSystem::ParameterCache paramCache;
        GetPropType<TypeTag, Properties::FluidState> fluidState;
        fluidState.setTemperature(this->spatialParams().temperatureAtPos(this->gridGeometry().bBoxMin()));
        fluidState.setPressure(FluidSystem::phase1Idx, initialPressure_);
        return FluidSystem::density(fluidState, paramCache, FluidSystem::phase1Idx);
    }

    std::shared_ptr<const FineLevelView> fineLevelView_;
    Scalar injectionRate_;
    Scalar initialPressure_;
};

} // end namespace Dumux

#endif
