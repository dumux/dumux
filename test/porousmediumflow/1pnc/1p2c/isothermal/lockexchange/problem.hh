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
 *        lighter (right) half at the same depth, released to develop a buoyancy-driven
 *        gravity current. Used to regression-test BoussinesqCVFEDarcyLaw against the
 *        full (variable-density) model, at a density variance small enough (a couple
 *        of percent) to be within the Boussinesq approximation's validity range.
 */

#ifndef DUMUX_LOCKEXCHANGE_TEST_PROBLEM_HH
#define DUMUX_LOCKEXCHANGE_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Lock-exchange test problem: closed (no-flow) domain, initial condition is a
 *        step function in tracer concentration across the vertical mid-plane (denser
 *        fluid on the left), hydrostatic within each half. Pressure is pinned at a
 *        single point (top-left corner, away from the initial interface) since the
 *        domain is otherwise fully closed (Neumann) for the pressure/mass equation.
 */
template <class TypeTag>
class LockExchangeTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    enum { pressureIdx = Indices::pressureIdx };

    //! The test is defined using mass fractions
    static_assert(!getPropValue<TypeTag, Properties::UseMoles>(), "This test uses mass fractions!");

public:
    LockExchangeTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        FluidSystem::init();
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    //! Closed box: no-flow everywhere, except a single pinned-pressure point (top-left
    //! corner, away from the initial interface at the domain's vertical mid-plane).
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if (globalPos[0] < eps_ && globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
            values.setDirichlet(pressureIdx);

        return values;
    }

    //! Pins the pressure to its (hydrostatic) initial value at the reference point.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Step-function initial condition: higher tracer concentration (denser
     *        fluid) in the left half of the domain, fresh (C=0) in the right half,
     *        each half hydrostatic with its own local density.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        const Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];
        const Scalar xMid = 0.5*(this->gridGeometry().bBoxMin()[0] + this->gridGeometry().bBoxMax()[0]);
        const Scalar concentration = (globalPos[0] < xMid) ? cLeft_ : 0.0;

        // hydrostatic within this column, using this test's exact (artificial) linear
        // density law -- always the "full" expression, independent of which
        // FluidSystem variant (full or Boussinesq) is actually active for this type
        // tag, so both variants start from the identical, physically consistent IC.
        const Scalar rho = referenceDensity_*(1.0 + concentration);

        PrimaryVariables priVars;
        priVars[pressureIdx] = pRef_ + rho*g_*depth;
        priVars[FluidSystem::soluteIdx] = concentration;
        return priVars;
    }

    // \}

    /*!
     * \brief The density deviation [kg/m^3] that drives buoyancy, for use with
     *        BoussinesqCVFEDarcyLaw (see dumux/flux/cvfe/boussinesqdarcyslaw.hh for
     *        the interface contract this method fulfills).
     *
     * Matches this test's exact (artificial) linear density law: rho0*C, so that
     * volVars.density() + buoyantDensity(volVars) reproduces the full model's
     * rho0*(1+C) exactly, whether or not volVars.density() itself already varies
     * with C (it doesn't, for the Boussinesq FluidSystem variant).
     */
    Scalar buoyantDensity(const VolumeVariables& volVars) const
    { return referenceDensity_*volVars.massFraction(0, FluidSystem::soluteIdx); }

private:
    static constexpr Scalar eps_ = 1e-6;
    static constexpr Scalar referenceDensity_ = 1000.0; // must match FluidSystem::referenceDensity()
    static constexpr Scalar cLeft_ = 0.01; // 1% density anomaly: within Boussinesq validity
    static constexpr Scalar g_ = 9.81;
    static constexpr Scalar pRef_ = 1e5;
};

} // end namespace Dumux

#endif
