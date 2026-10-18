// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Elder problem: dense fluid enters a closed porous layer through part of its top boundary
 */
#ifndef DUMUX_TEST_ELDER_PROBLEM_HH
#define DUMUX_TEST_ELDER_PROBLEM_HH

#include <cmath>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Right half of the Elder problem with the symmetry axis at the left boundary
 *
 * The solute mass fraction at the top is the smoothed step
 * \f$ X = \frac{1}{2} \left( 1 - \tanh\left(\frac{x - x_s}{\delta}\right) \right) \f$
 * with the half-width \f$x_s\f$ of the source and the transition length \f$\delta\f$, and the sharp step
 * for \f$\delta = 0\f$. The bottom has \f$X = 0\f$. No fluid crosses the boundaries.
 * The density deviation that causes buoyancy is \f$\Delta\varrho X\f$.
 */
template <class TypeTag>
class ElderProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GlobalPosition = typename GridGeometry::SubControlVolume::GlobalPosition;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    static constexpr int soluteIdx = FluidSystem::comp1Idx;
    static_assert(!getPropValue<TypeTag, Properties::UseMoles>(), "This test uses mass fractions");

public:
    ElderProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        FluidSystem::init();
        densityDifference_ = getParam<Scalar>("Problem.DensityDifference");
        sourceHalfWidth_ = getParam<Scalar>("Problem.SourceHalfWidth");
        transitionLength_ = getParam<Scalar>("Problem.SourceTransitionLength");
        height_ = gridGeometry->bBoxMax()[1] - gridGeometry->bBoxMin()[1];
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (onTop_(globalPos) || onBottom_(globalPos))
            values.setDirichlet(soluteIdx);
        if (onTop_(globalPos) && globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
            values.setDirichlet(Indices::pressureIdx);
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }

    //! Hydrostatic pressure of the fluid without solute, the source concentration at the top and no solute elsewhere
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        const Scalar z = globalPos[1] - this->gridGeometry().bBoxMin()[1];
        const Scalar rho = FluidSystem::density(this->spatialParams().temperatureAtPos(globalPos), 1e5);
        const Scalar g = -this->spatialParams().gravity(globalPos)[1];

        PrimaryVariables values;
        values[Indices::pressureIdx] = 1e5 + rho*g*(height_ - z);
        values[soluteIdx] = onTop_(globalPos) ? sourceMassFraction(globalPos[0]) : 0.0;
        return values;
    }

    //! The density deviation that causes buoyancy
    Scalar buoyantDensity(const VolumeVariables& volVars, int phaseIdx) const
    { return densityDifference_*volVars.massFraction(phaseIdx, soluteIdx); }

    //! The density deviation at solute mass fraction 1
    Scalar densityDifference() const
    { return densityDifference_; }

    //! The solute mass fraction at the top boundary as a function of the horizontal coordinate
    Scalar sourceMassFraction(Scalar x) const
    {
        const Scalar distance = x - this->gridGeometry().bBoxMin()[0] - sourceHalfWidth_;
        if (transitionLength_ > 0.0)
            return 0.5*(1.0 - std::tanh(distance/transitionLength_));
        return std::abs(distance) < eps_ ? 0.5 : (distance < 0.0 ? 1.0 : 0.0);
    }

    bool onTop(const GlobalPosition& globalPos) const
    { return onTop_(globalPos); }

    bool onBottom(const GlobalPosition& globalPos) const
    { return onBottom_(globalPos); }

private:
    bool onTop_(const GlobalPosition& globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    bool onBottom_(const GlobalPosition& globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    static constexpr Scalar eps_ = 1e-6;
    Scalar densityDifference_, sourceHalfWidth_, transitionLength_, height_;
};

} // end namespace Dumux

#endif
