// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Horton-Rogers-Lapwood problem: onset of convection in a porous layer with dense fluid on top
 */
#ifndef DUMUX_TEST_HORTON_ROGERS_LAPWOOD_PROBLEM_HH
#define DUMUX_TEST_HORTON_ROGERS_LAPWOOD_PROBLEM_HH

#include <cmath>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Closed box with solute mass fraction 1 at the top and 0 at the bottom
 *
 * The initial condition is the linear profile at rest plus the normal mode
 * \f$ A \cos(k x) \sin(\pi z / H) \f$ with \f$ k = \pi / L \f$.
 * The density deviation that causes buoyancy is \f$\Delta\varrho X\f$.
 */
template <class TypeTag>
class HortonRogersLapwoodProblem : public PorousMediumFlowProblem<TypeTag>
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
    HortonRogersLapwoodProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        FluidSystem::init();
        densityDifference_ = getParam<Scalar>("Problem.DensityDifference");
        amplitude_ = getParam<Scalar>("Problem.Amplitude");
        length_ = gridGeometry->bBoxMax()[0] - gridGeometry->bBoxMin()[0];
        height_ = gridGeometry->bBoxMax()[1] - gridGeometry->bBoxMin()[1];
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (onTop_(globalPos) || onBottom_(globalPos))
            values.setDirichlet(soluteIdx);
        if (onBottom_(globalPos) && globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
            values.setDirichlet(Indices::pressureIdx);
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        auto values = initialAtPos(globalPos);
        values[soluteIdx] = onTop_(globalPos) ? 1.0 : 0.0;
        return values;
    }

    //! Linear profile with hydrostatic pressure plus the normal mode
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        const Scalar z = globalPos[1] - this->gridGeometry().bBoxMin()[1];
        const Scalar rho = FluidSystem::density(this->spatialParams().temperatureAtPos(globalPos), 1e5);
        const Scalar g = -this->spatialParams().gravity(globalPos)[1];

        PrimaryVariables values;
        values[Indices::pressureIdx] = 1e5 + g*(rho*(height_ - z) + 0.5*densityDifference_*(height_*height_ - z*z)/height_);
        values[soluteIdx] = baseProfile(globalPos) + amplitude_*mode(globalPos);
        return values;
    }

    //! The density deviation that causes buoyancy
    Scalar buoyantDensity(const VolumeVariables& volVars, int phaseIdx) const
    { return densityDifference_*volVars.massFraction(phaseIdx, soluteIdx); }

    //! The density deviation at solute mass fraction 1
    Scalar densityDifference() const
    { return densityDifference_; }

    Scalar baseProfile(const GlobalPosition& globalPos) const
    { return (globalPos[1] - this->gridGeometry().bBoxMin()[1])/height_; }

    Scalar mode(const GlobalPosition& globalPos) const
    {
        using std::cos; using std::sin;
        const Scalar x = globalPos[0] - this->gridGeometry().bBoxMin()[0];
        const Scalar z = globalPos[1] - this->gridGeometry().bBoxMin()[1];
        return cos(M_PI*x/length_)*sin(M_PI*z/height_);
    }

private:
    bool onTop_(const GlobalPosition& globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    bool onBottom_(const GlobalPosition& globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    static constexpr Scalar eps_ = 1e-8;
    Scalar densityDifference_, amplitude_, length_, height_;
};

} // end namespace Dumux

#endif
