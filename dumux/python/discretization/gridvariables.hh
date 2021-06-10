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

#ifndef DUMUX_PYTHON_DISCRETIZATION_GRIDVARIABLES_HH
#define DUMUX_PYTHON_DISCRETIZATION_GRIDVARIABLES_HH

#include <memory>
#include <dune/common/classname.hh>
#include <dune/istl/bvector.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/common/typeregistry.hh>

namespace Dumux::Python {


// see python/dumux/discretization/__init__.py for how this is used for JIT compilation
template <class GV, class... Options>
void registerGridVariables(pybind11::handle scope, pybind11::class_<GV, Options...> cls)
{
    using pybind11::operator""_a;

    using Problem = typename GV::GridVolumeVariables::Problem;
    using PrimaryVariables = typename GV::GridVolumeVariables::VolumeVariables::PrimaryVariables;
    using GridGeometry = typename GV::GridGeometry;
    using SolutionVector = Dune::BlockVector<PrimaryVariables>;

    cls.def(pybind11::init([](std::shared_ptr<const Problem> problem,
                              std::shared_ptr<const GridGeometry> gridGeometry){
        return std::make_shared<GV>(problem, gridGeometry);
    }));

    // cls.def("update", &GG::update);
    cls.def("curGridVolVars", [](GV& self) { return self.curGridVolVars(); });

    cls.def("init", [](GV& self, const SolutionVector& sol) { return self.init(sol); });

//     cls.def("numDofs", &GG::numDofs);
//     cls.def("numScv", &GG::numScv);
//     cls.def("numScvf", &GG::numScvf);

}

} // namespace Dumux::Python

#endif
