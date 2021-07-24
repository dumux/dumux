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
#include <dune/istl/bvector.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/common/typeregistry.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux::Python {

namespace Impl {

template<class ElementSolution>
void registerElementSolution(pybind11::handle scope)
{
    using namespace Dune::Python;

    auto [cls, addedToRegistry] = insertClass<ElementSolution>(
        scope, "ElementSolution",
        GenerateTypeName(Dune::className<ElementSolution>()),
        IncludeFiles{"dumux/discretization/elementsolution.hh"}
    );

    if (addedToRegistry)
    {
        using pybind11::operator""_a;

        cls.def("__getitem__", [](const ElementSolution& self, std::size_t i){
            if (i >= self.size())
                throw pybind11::index_error();
            return self[i];
        });

        cls.def_property_readonly("size", &ElementSolution::size);
    }
}

} // end namespace Impl


// see python/dumux/discretization/__init__.py for how this is used for JIT compilation
template <class GV, class... Options>
void registerGridVariables(pybind11::handle scope, pybind11::class_<GV, Options...> cls)
{
    using pybind11::operator""_a;

    using Problem = typename GV::GridVolumeVariables::Problem;
    using PrimaryVariables = typename GV::GridVolumeVariables::VolumeVariables::PrimaryVariables;
    using GridGeometry = typename GV::GridGeometry;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SolutionVector = Dune::BlockVector<PrimaryVariables>;

    cls.def(pybind11::init([](std::shared_ptr<const Problem> problem,
                              std::shared_ptr<const GridGeometry> gridGeometry){
        return std::make_shared<GV>(problem, gridGeometry);
    }));

    cls.def("init", [](GV& self, const SolutionVector& sol) { return self.init(sol); });
    cls.def("advanceTimeStep", &GV::advanceTimeStep);
    cls.def_property_readonly("curGridVolVars", [](GV& self) { return self.curGridVolVars(); });
    cls.def_property_readonly("gridFluxVarsCache", [](GV& self) { return self.gridFluxVarsCache(); });
    cls.def_property_readonly("prevGridVolVars", [](GV& self) { return self.prevGridVolVars(); });
    cls.def_property_readonly("gridGeometry", &GV::gridGeometry);

    cls.def("updateAfterGridAdaption", [](GV& self, const SolutionVector& sol){
        return self.updateAfterGridAdaption(sol);
    });

    cls.def("resetTimeStep", [](GV& self, const SolutionVector& sol){
        return self.resetTimeStep(sol);
    });

    cls.def("update", [](GV& self, const SolutionVector& sol, const bool forceFluxCacheUpdate = false){
        return self.update(sol, forceFluxCacheUpdate);
    });

    using ElementSolution = std::decay_t<decltype(elementSolution(std::declval<Element>(),
                                                                  std::declval<SolutionVector>(),
                                                                  std::declval<GridGeometry>()))>;
    Impl::registerElementSolution<ElementSolution>(scope);
}

} // namespace Dumux::Python

#endif
