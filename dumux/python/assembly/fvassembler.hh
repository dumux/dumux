// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief TODO: docme!
 */

#ifndef DUMUX_PYTHON_COMMON_FVASSEMBLER_HH
#define DUMUX_PYTHON_COMMON_FVASSEMBLER_HH

#include <limits>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <dumux/common/timeloop.hh>

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
    using Scalar = typename FVAssembler::Scalar;

    static_assert(std::is_same_v<GridGeometry, typename Problem::GridGeometry>);
    cls.def(pybind11::init([](std::shared_ptr<const Problem> problem,
                              std::shared_ptr<const GridGeometry> gridGeometry,
                              std::shared_ptr<GridVariables> gridVariables){
        return std::make_shared<FVAssembler>(problem, gridGeometry, gridVariables);
    }));

    // instationary assembly; the time loop is created here so that none has to be passed
    // across the language boundary
    cls.def("assembleJacobianAndResidual",
        [](FVAssembler& self, const SolutionVector& curSol,
           const SolutionVector& prevSol, Scalar dt){
            auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(
                Scalar(0), dt, std::numeric_limits<Scalar>::max(), /*verbose=*/false);
            self.setTimeLoop(timeLoop);
            self.setPreviousSolution(prevSol);
            self.assembleJacobianAndResidual(curSol);
        }, "curSol"_a, "prevSol"_a, "dt"_a);

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
