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

#ifndef DUMUX_PYTHON_POROUSMEDIUMFLOW_PROBLEM_HH
#define DUMUX_PYTHON_POROUSMEDIUMFLOW_PROBLEM_HH


#include <dune/python/pybind11/pybind11.h>

#include <dumux/python/common/fvproblem.hh>

namespace Dumux::Python {

/*!
 * \ingroup PorousmediumflowModels
 * \brief A C++ wrapper for a Python PorousMediumFlow problem
 */
template<class GridGeometry_, class SpatialParams_, class PrimaryVariables, bool enableInternalDirichletConstraints>
class PorousMediumFlowProblem
: public Dumux::Python::FVProblem<GridGeometry_, SpatialParams_, PrimaryVariables, enableInternalDirichletConstraints>
{
    using ParentType = Dumux::Python::FVProblem<GridGeometry_, SpatialParams_, PrimaryVariables, enableInternalDirichletConstraints>;
public:
    using GridGeometry = GridGeometry_;
    using SpatialParams = SpatialParams_;
    using Scalar = typename PrimaryVariables::value_type;
    using NumEqVector = Dune::FieldVector<Scalar, PrimaryVariables::dimension>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethods::box;
    static constexpr std::size_t numEq = static_cast<std::size_t>(PrimaryVariables::dimension);
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::dimension>;

    PorousMediumFlowProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                            std::shared_ptr<const SpatialParams> spatialParams,
                            pybind11::object pyProblem)
    : ParentType(gridGeometry, spatialParams, pyProblem)
    , spatialParams_(spatialParams)
    {}

    const SpatialParams& spatialParams() const
    { return *spatialParams_; }

private:
    std::shared_ptr<const SpatialParams> spatialParams_;
};

// Python wrapper for the above FVProblem C++ class
template<class Problem, class... options>
void registerPorousMediumFlowProblem(pybind11::handle scope, pybind11::class_<Problem, options...> cls)
{
    using pybind11::operator""_a;
    using namespace Dune::Python;

    using GridGeometry = typename Problem::GridGeometry;
    using SpatialParams = typename Problem::SpatialParams;
    cls.def(pybind11::init([](std::shared_ptr<const GridGeometry> gridGeometry,
                              std::shared_ptr<SpatialParams> spatialParams,
                              pybind11::object p){
        return std::make_shared<Problem>(gridGeometry, spatialParams, p);
    }));

    cls.def_property_readonly("name", &Problem::name);
    cls.def_property_readonly("numEq", [](Problem&){ return Problem::numEq; });

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

    if constexpr (Problem::isBox)
    {
        using SCV = typename Problem::SubControlVolume;
        cls.def("boundaryTypes", pybind11::overload_cast<const Element&, const SCV&>(&Problem::boundaryTypes, pybind11::const_), "element"_a, "scv"_a);
        cls.def("dirichlet", pybind11::overload_cast<const Element&, const SCV&>(&Problem::dirichlet, pybind11::const_), "element"_a, "scv"_a);
    }
    else
    {
        using SCVF = typename Problem::SubControlVolumeFace;
        cls.def("boundaryTypes", pybind11::overload_cast<const Element&, const SCVF&>(&Problem::boundaryTypes, pybind11::const_), "element"_a, "scvf"_a);
        cls.def("dirichlet", pybind11::overload_cast<const Element&, const SCVF&>(&Problem::dirichlet, pybind11::const_), "element"_a, "scvf"_a);
    }

    cls.def("neumann", &Problem::template neumann<decltype(std::ignore), decltype(std::ignore)>);
    cls.def("source", &Problem::template source<decltype(std::ignore)>);
    cls.def("sourceAtPos", &Problem::sourceAtPos);
    cls.def("initial", &Problem::template initial<Element>);
    cls.def("initial", &Problem::template initial<Vertex>);
    cls.def("gridGeometry", &Problem::gridGeometry);
    cls.def("spatialParams", &Problem::spatialParams);
}

} // end namespace Dumux::Python

#endif
