// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief The system u0' = -u0 + u1, 0 = sin(t) - u1 of the time stepping tests, whose second
 *        equation has no storage term
 */
#ifndef DUMUX_TEST_TIMESTEPPING_CONSTRAINED_ODE_HH
#define DUMUX_TEST_TIMESTEPPING_CONSTRAINED_ODE_HH

#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/experimental/ode/odesolver.hh>

namespace Dumux::Test {

//! u0' = -u0 + u1 with the constraint 0 = sin(t) - u1 (no storage term in the second equation)
class ConstrainedODE
{
public:
    using Scalar = double;
    using SolutionVector = Dune::FieldVector<Scalar, 2>;
    using ResidualType = SolutionVector;
    using JacobianMatrix = Dune::FieldMatrix<Scalar, 2, 2>;
    using Variables = Experimental::ODEVariables<SolutionVector>;

    void accumulation(const Variables& vars, ResidualType& value) const
    {
        value[0] = vars.dofs()[0];
        value[1] = 0.0;
    }

    void accumulationJacobian(const Variables&, JacobianMatrix& jacobian) const
    {
        jacobian = 0.0;
        jacobian[0][0] = 1.0;
    }

    void rhs(const Variables& vars, ResidualType& rhs) const
    {
        using std::sin;
        rhs[0] = -vars.dofs()[0] + vars.dofs()[1];
        rhs[1] = sin(vars.independentVariableLevel().current()) - vars.dofs()[1];
    }

    void rhsJacobian(const Variables&, JacobianMatrix& jacobian) const
    {
        jacobian = 0.0;
        jacobian[0][0] = -1.0;
        jacobian[0][1] = 1.0;
        jacobian[1][1] = -1.0;
    }
};

//! u0(t) of the constrained system for u0(0) = 0
inline double constrainedExact(const double t)
{
    using std::sin; using std::cos; using std::exp;
    return 0.5*(sin(t) - cos(t) + exp(-t));
}

} // end namespace Dumux::Test

#endif
