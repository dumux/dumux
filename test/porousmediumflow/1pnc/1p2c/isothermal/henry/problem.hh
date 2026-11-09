// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief The Henry problem of Fahs et al. \cite Fahs2016
 */
#ifndef DUMUX_TEST_HENRY_FAHS_PROBLEM_HH
#define DUMUX_TEST_HENRY_FAHS_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Seawater intrusion into a confined aquifer (Henry \cite Henry1964) as defined by Fahs et al. \cite Fahs2016
 */
template <class TypeTag>
class HenryFahsProblem : public PorousMediumFlowProblem<TypeTag>
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
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    static constexpr int soluteIdx = FluidSystem::comp1Idx;
    static_assert(!getPropValue<TypeTag, Properties::UseMoles>(), "This test uses mass fractions");

public:
    HenryFahsProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        FluidSystem::init();
        densityDifference_ = getParam<Scalar>("Problem.DensityDifference");
        length_ = gridGeometry->bBoxMax()[0] - gridGeometry->bBoxMin()[0];
        depth_ = gridGeometry->bBoxMax()[1] - gridGeometry->bBoxMin()[1];

        const auto& pos = gridGeometry->bBoxMin();
        const Scalar temperature = this->spatialParams().temperatureAtPos(pos);
        density_ = FluidSystem::density(temperature, referencePressure_);
        gravity_ = -this->spatialParams().gravity(pos)[1];

        // uniform freshwater inflow (mass flux, depth is the length of the inflow boundary)
        inflowFlux_ = density_*getParam<Scalar>("Problem.FreshwaterInflow")/depth_;
    }

    //! Sea: seawater at rest (Dirichlet), elsewhere Neumann
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (onRight_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }

    //! Uniform inflow of freshwater (negative: into the domain) on the left, no flow elsewhere
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);
        if (onLeft_(globalPos))
            values[FluidSystem::comp0Idx] = -inflowFlux_;
        return values;
    }

    //! Seawater at rest
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = referencePressure_ + (density_ + densityDifference_)*gravity_*depthBelowTop_(globalPos);
        values[soluteIdx] = 1.0;
        return values;
    }

    //! The density deviation that causes buoyancy
    Scalar buoyantDensity(const VolumeVariables& volVars, int phaseIdx) const
    { return densityDifference_*volVars.massFraction(phaseIdx, soluteIdx); }

private:
    Scalar depthBelowTop_(const GlobalPosition& globalPos) const
    { return this->gridGeometry().bBoxMax()[1] - globalPos[1]; }

    bool onLeft_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_*length_; }

    bool onRight_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_*length_; }

    static constexpr Scalar eps_ = 1e-8;
    static constexpr Scalar referencePressure_ = 1e5;
    Scalar densityDifference_, length_, depth_, density_, gravity_, inflowFlux_;
};

} // end namespace Dumux

#endif
