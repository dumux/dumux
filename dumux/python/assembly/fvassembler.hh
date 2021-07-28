// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \brief TODO: docme!
 */

#ifndef DUMUX_PYTHON_COMMON_FVASSEMBLER_HH
#define DUMUX_PYTHON_COMMON_FVASSEMBLER_HH

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

namespace Dumux::Python {

// Python wrapper for the FVAssembler C++ class
template<class FVAssembler, class... options>
void registerFVAssembler(pybind11::handle scope, pybind11::class_<FVAssembler, options...> cls)
{
    using pybind11::operator""_a;

    using Problem = typename FVAssembler::Problem;
    using GridGeometry = typename FVAssembler::GridGeometry;
    using GridVariables = typename FVAssembler::GridVariables;
    using SolutionVector = typename FVAssembler::ResidualType;

    static_assert(std::is_same_v<GridGeometry, typename Problem::GridGeometry>);
    cls.def(pybind11::init([](std::shared_ptr<const Problem> problem,
                              std::shared_ptr<const GridGeometry> gridGeometry,
                              std::shared_ptr<GridVariables> gridVariables){
        return std::make_shared<FVAssembler>(problem, gridGeometry, gridVariables);
    }));

    // TODO assembler with time loop

    cls.def_property_readonly("numDofs", &FVAssembler::numDofs);
    cls.def_property_readonly("problem", &FVAssembler::problem);
    cls.def_property_readonly("gridGeometry", &FVAssembler::gridGeometry);
    cls.def_property_readonly("gridView", &FVAssembler::gridView);
    cls.def_property_readonly("jacobian", &FVAssembler::jacobian);
    cls.def_property_readonly("residual", &FVAssembler::residual);
    cls.def_property_readonly("prevSol", &FVAssembler::prevSol);
    cls.def_property_readonly("isStationaryProblem", &FVAssembler::isStationaryProblem);
    cls.def_property_readonly("gridVariables", [](FVAssembler& self) { return self.gridVariables(); });

    cls.def("updateGridVariables", [](FVAssembler& self, const SolutionVector& curSol){
        self.updateGridVariables(curSol);
    });

    cls.def("assembleResidual", [](FVAssembler& self, const SolutionVector& curSol){
        self.assembleResidual(curSol);
    });

    cls.def("assembleJacobianAndResidual", [](FVAssembler& self, const SolutionVector& curSol){
        self.assembleJacobianAndResidual(curSol);
    });
}

} // end namespace Dumux::Python

#endif
