// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Lock-exchange test: a closed box, initially divided into a denser (left) and
 *        lighter (right) half, develops a buoyancy-driven gravity current.
 */
#ifndef DUMUX_LOCKEXCHANGE_TEST_PROBLEM_HH
#define DUMUX_LOCKEXCHANGE_TEST_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Lock-exchange test problem with the Boussinesq approximation
 *
 * The domain is closed. The solute mass fraction is a step across the vertical mid-plane
 * (solute on the left), and the density deviation is \f$\Delta\varrho = \varrho X\f$.
 * The pressure is fixed at the top-left corner.
 */
template <class TypeTag>
class LockExchangeTestProblem : public PorousMediumFlowProblem<TypeTag>
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
    static_assert(!getPropValue<TypeTag, Properties::UseMoles>(), "This test uses mass fractions!");

public:
    LockExchangeTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        FluidSystem::init();
        massFractionLeft_ = getParam<Scalar>("Problem.MassFractionLeft");
        g_ = -this->spatialParams().gravity(GlobalPosition(0.0))[1];
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (globalPos[0] < eps_ && globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
            values.setDirichlet(Indices::pressureIdx);
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }

    //! Hydrostatic pressure within each half
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        const Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];
        const Scalar xMid = 0.5*(this->gridGeometry().bBoxMin()[0] + this->gridGeometry().bBoxMax()[0]);
        const Scalar massFraction = (globalPos[0] < xMid) ? massFractionLeft_ : 0.0;

        PrimaryVariables priVars;
        priVars[Indices::pressureIdx] = 1e5 + referenceDensity_()*(1.0 + massFraction)*g_*depth;
        priVars[soluteIdx] = massFraction;
        return priVars;
    }

    //! The density deviation that causes buoyancy
    Scalar buoyantDensity(const VolumeVariables& volVars, int phaseIdx) const
    { return volVars.density(phaseIdx)*volVars.massFraction(phaseIdx, soluteIdx); }

private:
    Scalar referenceDensity_() const
    { return FluidSystem::density(this->spatialParams().temperatureAtPos(GlobalPosition(0.0)), 1e5); }

    static constexpr Scalar eps_ = 1e-6;
    Scalar massFractionLeft_;
    Scalar g_;
};

} // end namespace Dumux

#endif
