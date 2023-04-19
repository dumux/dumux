// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
    using SolutionVector = typename FVAssembler::SolutionVector;

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

    cls.def("assembleResidual", [](FVAssembler& self, const SolutionVector& curSol){
        self.assembleResidual(curSol);
    });

    cls.def("assembleJacobianAndResidual", [](FVAssembler& self, const SolutionVector& curSol){
        self.assembleJacobianAndResidual(curSol);
    });
}

} // end namespace Dumux::Python

#endif
