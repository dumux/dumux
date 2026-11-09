// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief One-sided convection: a porous layer below a boundary that is saturated with a dense solute
 */
#ifndef DUMUX_TEST_ONE_SIDED_CONVECTION_PROBLEM_HH
#define DUMUX_TEST_ONE_SIDED_CONVECTION_PROBLEM_HH

#include <cmath>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Closed box with solute mass fraction 1 at the top and no solute flux through the other boundaries
 *
 * With depth \f$ d \f$ below the top, dimensionless time \f$ t \f$ and the Rayleigh number \f$ Ra \f$, the initial
 * condition at \f$ t_0 \f$ is the diffusive profile \f$ \mathrm{erfc}(\zeta) \f$, \f$ \zeta = d / (2 H \sqrt{t_0/Ra}) \f$,
 * plus the mode \f$ A \cos(k x) \zeta \exp(-\zeta^2) \f$ with \f$ k = \pi / L \f$.
 * The density deviation that causes buoyancy is \f$\Delta\varrho X\f$.
 */
template <class TypeTag>
class OneSidedConvectionProblem : public PorousMediumFlowProblem<TypeTag>
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
    OneSidedConvectionProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        FluidSystem::init();
        densityDifference_ = getParam<Scalar>("Problem.DensityDifference");
        amplitude_ = getParam<Scalar>("Problem.Amplitude");
        initialTime_ = getParam<Scalar>("Problem.InitialTime");
        length_ = gridGeometry->bBoxMax()[0] - gridGeometry->bBoxMin()[0];
        height_ = gridGeometry->bBoxMax()[1] - gridGeometry->bBoxMin()[1];
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (onTop_(globalPos))
            values.setDirichlet(soluteIdx);
        if (onBottom_(globalPos) && globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
            values.setDirichlet(Indices::pressureIdx);
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        auto values = initialAtPos(globalPos);
        values[soluteIdx] = 1.0;
        return values;
    }

    //! Diffusive profile at the initial time with hydrostatic pressure, plus the mode
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        using std::erfc; using std::exp; using std::sqrt;
        const Scalar d = depth(globalPos);
        const Scalar b = 2.0*height_*sqrt(initialTime_/rayleigh_);
        const Scalar rho = FluidSystem::density(this->spatialParams().temperatureAtPos(globalPos), 1e5);
        const Scalar g = -this->spatialParams().gravity(globalPos)[1];
        const Scalar integratedProfile = d*erfc(d/b) + b/sqrt(M_PI)*(1.0 - exp(-d*d/(b*b)));

        PrimaryVariables values;
        values[Indices::pressureIdx] = 1e5 + g*(rho*d + densityDifference_*integratedProfile);
        values[soluteIdx] = baseProfile(globalPos, initialTime_) + amplitude_*mode(globalPos);
        return values;
    }

    //! The density deviation that causes buoyancy
    Scalar buoyantDensity(const VolumeVariables& volVars, int phaseIdx) const
    { return densityDifference_*volVars.massFraction(phaseIdx, soluteIdx); }

    //! The density deviation at solute mass fraction 1
    Scalar densityDifference() const
    { return densityDifference_; }

    //! The dimensionless time of the initial condition
    Scalar initialTime() const
    { return initialTime_; }

    //! Set the Rayleigh number used by the initial condition
    void setRayleighNumber(Scalar rayleigh)
    { rayleigh_ = rayleigh; }

    Scalar depth(const GlobalPosition& globalPos) const
    { return this->gridGeometry().bBoxMax()[1] - globalPos[1]; }

    //! Diffusive profile of the semi-infinite layer at dimensionless time t
    Scalar baseProfile(const GlobalPosition& globalPos, Scalar t) const
    {
        using std::erfc; using std::sqrt;
        return erfc(depth(globalPos)/(2.0*height_*sqrt(t/rayleigh_)));
    }

    Scalar mode(const GlobalPosition& globalPos) const
    {
        using std::cos; using std::exp; using std::sqrt;
        const Scalar x = globalPos[0] - this->gridGeometry().bBoxMin()[0];
        const Scalar zeta = depth(globalPos)/(2.0*height_*sqrt(initialTime_/rayleigh_));
        return cos(M_PI*x/length_)*zeta*exp(-zeta*zeta);
    }

private:
    bool onTop_(const GlobalPosition& globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    bool onBottom_(const GlobalPosition& globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    static constexpr Scalar eps_ = 1e-8;
    Scalar densityDifference_, amplitude_, initialTime_, length_, height_;
    Scalar rayleigh_ = 0.0;
};

} // end namespace Dumux

#endif
