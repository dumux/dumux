// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePNCTests
 * \brief This directory implements the Henry problem benchmark of Fahs et al. (2016,
 *        WRR, doi:10.1002/2016WR019288), "Test Case 1": the classic (purely diffusive)
 *        Henry (1964) configuration. Used to validate this model's flow/transport
 *        implementation against a published semianalytical solution (see
 *        validate_fahs2016.py and Appendix D, Table D1 of the paper).
 *
 * The sea (right) boundary is the *original* Henry (1964)/Fahs et al. (2016)
 * condition: a simple Dirichlet condition on both pressure (hydrostatic) and
 * concentration (fixed at seawater salinity, their "c=1"), not the Segol/Voss mixed
 * inflow/outflow condition sometimes used for this problem elsewhere in the
 * literature -- deliberately the simpler boundary condition, matching the reference
 * solution being validated against.
 */

#ifndef DUMUX_HENRY_FAHS_TEST_PROBLEM_HH
#define DUMUX_HENRY_FAHS_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief The classic (purely diffusive) Henry problem (Henry 1964; Fahs et al. 2016
 *        Test Case 1).
 */
template <class TypeTag>
class HenryFahsTestProblem : public PorousMediumFlowProblem<TypeTag>
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
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    // copy pressure index for convenience
    enum { pressureIdx = Indices::pressureIdx };

    //! The test is defined using mass fractions
    static_assert(!getPropValue<TypeTag, Properties::UseMoles>(), "This test uses mass fractions!");

public:
    HenryFahsTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        FluidSystem::init();

        domainHeight_ = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];

        // freshwater inflow rate [m^3/d], divided by the (2D, unit-thickness) left
        // boundary's area to get a specific discharge [m/s]
        const Scalar inflowRate = getParam<Scalar>("Problem.InflowRate"); // [m^3/d]
        freshwaterSpecificDischarge_ = inflowRate/86400.0/domainHeight_;  // [m/s]
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Left: Neumann everywhere (specified freshwater inflow, see
     *        neumannAtPos()). Right: Dirichlet everywhere -- hydrostatic pressure and
     *        fixed (seawater) concentration, the original Henry (1964)/Fahs et al.
     *        (2016) sea boundary condition.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
            values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Hydrostatic pressure (seawater reference density) and fixed seawater
     *        concentration at the right boundary.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        const Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];
        PrimaryVariables values;
        values[pressureIdx] = pRef_ + FluidSystem::rhoSeawater()*g_*depth;
        values[FluidSystem::soluteIdx] = FluidSystem::seawaterMassFraction();
        return values;
    }

    /*!
     * \brief Specified freshwater inflow (fixed rate, zero salt) at the left boundary;
     *        no-flow everywhere else (the right boundary is Dirichlet, so this is
     *        only ever evaluated on the left and on the impermeable top/bottom).
     */
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);

        if (globalPos[0] < eps_)
        {
            // negative sign as this is an influx; no salt is carried in
            values[FluidSystem::solventIdx] = -freshwaterSpecificDischarge_*FluidSystem::rhoFresh();
            values[FluidSystem::soluteIdx] = 0.0;
        }

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Initial condition: the domain is initially filled with seawater,
     *        hydrostatic (using the seawater reference density).
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        const Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];
        PrimaryVariables values;
        values[pressureIdx] = pRef_ + FluidSystem::rhoSeawater()*g_*depth;
        values[FluidSystem::soluteIdx] = FluidSystem::seawaterMassFraction();
        return values;
    }

    // \}

private:
    static constexpr Scalar eps_ = 1e-6;
    static constexpr Scalar g_ = 9.81;
    static constexpr Scalar pRef_ = 1e5;

    Scalar domainHeight_;
    Scalar freshwaterSpecificDischarge_;
};

} // end namespace Dumux

#endif
