// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief The problem for the analytic-Jacobian assembly test
 */
#ifndef DUMUX_TEST_EXPERIMENTAL_ASSEMBLY_PROBLEM_HH
#define DUMUX_TEST_EXPERIMENTAL_ASSEMBLY_PROBLEM_HH

#include <cmath>

#include <dumux/common/problem.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \brief Nonlinear diffusion on the unit square, no-flux everywhere
 *
 * The boundary is left at zero flux and no constraints are set: the test compares the assembled
 * Jacobian, and every row should be a balance rather than a constraint the assembler writes
 * identically whichever differentiation method produced it.
 */
template<class TypeTag>
class NonlinearDiffusionTestProblem : public Experimental::Problem<TypeTag>
{
    using ParentType = Experimental::Problem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using NumEqVector = typename ParentType::Traits::NumEqVector;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<1>;

    using ParentType::ParentType;

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes values;
        values.setAllFluxBoundary();
        return values;
    }

    //! zero flux, so the boundary integral has no derivative to contribute
    template<class ElementVariables, class IpData>
    NumEqVector boundaryFlux(const ElementDiscretization&,
                             const ElementVariables&,
                             const IpData&) const
    { return NumEqVector(0.0); }

    //! a position-dependent source, so the interior residual is not identically zero
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        using std::sin;
        return NumEqVector(sin(2.0*M_PI*globalPos[0])*sin(2.0*M_PI*globalPos[1]));
    }
};

} // end namespace Dumux

#endif
