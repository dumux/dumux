// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Hydrostatic test problem for the one-phase model
 */

#ifndef DUMUX_TEST_ONEP_HYDROSTATIC_PROBLEM_HH
#define DUMUX_TEST_ONEP_HYDROSTATIC_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Hydrostatic test problem for the one-phase model.
 *
 * Gravity is enabled, a hydrostatic pressure is prescribed on the left boundary and all
 * other boundaries are impermeable. The exact solution is the hydrostatic pressure
 * \f$ p = p_0 + \varrho \mathbf{g} \cdot \mathbf{x} \f$ with zero flux everywhere,
 * independent of the grid and the permeability distribution.
 */
template<class TypeTag>
class OnePHydrostaticProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;

public:
    OnePHydrostaticProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(exactPressure(globalPos)); }

    //! The exact (hydrostatic) pressure
    Scalar exactPressure(const GlobalPosition& globalPos) const
    {
        const Scalar rho = FluidSystem::density(this->spatialParams().temperatureAtPos(globalPos), p0_);
        return p0_ + rho*(this->spatialParams().gravity(globalPos)*globalPos);
    }

private:
    static constexpr Scalar p0_ = 1.0e5;
    static constexpr Scalar eps_ = 1.0e-6;
};

} // end namespace Dumux

#endif
