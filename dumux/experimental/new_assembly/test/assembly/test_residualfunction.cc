// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Tests for the `ResidualFunction` class.
 */
#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/initialize.hh>

#include <dumux/experimental/new_assembly/dumux/assembly/residualfunction.hh>

using Scalar = double;
using Residual = Dune::FieldVector<Scalar, 1>;
using Jacobian = Dune::FieldMatrix<Scalar, 1, 1>;
using Derivative = Jacobian;

// assembles the equation x^2 - 5 = 0
struct MyAssembler
{
    using Variables = Scalar;

    void assembleResidual(Residual& r, const Variables& v) const
    { r[0] = v*v - 5.0; }

    void assembleJacobianAndResidual(Jacobian& j, Residual& r, const Variables& v) const
    {
        assembleResidual(r, v);
        j[0][0] = 2.0*v;
    }
};

int main (int argc, char *argv[])
{
    Dumux::initialize(argc, argv);

    Dumux::ResidualFunction residual{
        std::make_shared<const MyAssembler>(),
        std::make_shared<Jacobian>(),
        std::make_shared<Residual>()
    };

    if (std::abs(residual.evaluateAt(0.0) + 5.0) > 1e-6)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected residual");
    if (std::abs(residual.evaluateAt(2.0) + 1.0) > 1e-6)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected residual");

    const auto linearization = residual.linearizeAt(2.0);
    if (std::abs(linearization.value()[0] + 1.0) > 1e-6)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected residual");
    if (std::abs(linearization.derivative()[0][0] - 4.0) > 1e-6)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected derivative");

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
