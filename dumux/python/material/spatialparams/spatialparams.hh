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

#ifndef DUMUX_PYTHON_MATERIAL_SPATIAL_PARAMS_HH
#define DUMUX_PYTHON_MATERIAL_SPATIAL_PARAMS_HH

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <dumux/porousmediumflow/fvspatialparams.hh>

namespace Dumux::Python {

template<class GridGeometry, class Scalar, class PT>
class FVSpatialParamsOneP
: public Dumux::FVSpatialParamsOneP<GridGeometry, Scalar, FVSpatialParamsOneP<GridGeometry, Scalar, PT>>
{
    using ThisType = FVSpatialParamsOneP<GridGeometry, Scalar, PT>;
    using ParentType = Dumux::FVSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;

public:

    using PermeabilityType = PT;

    FVSpatialParamsOneP(std::shared_ptr<const GridGeometry> gridGeometry,
                        pybind11::object pySpatialParams)
    : ParentType(gridGeometry)
    , pySpatialParams_(pySpatialParams)
    {}

    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        if (pybind11::hasattr(pySpatialParams_, "permeability"))
            return pySpatialParams_.attr("permeability")(element, scv, elemSol).template cast<PermeabilityType>();
        else
            return pySpatialParams_.attr("permeabilityAtPos")(scv.center()).template cast<PermeabilityType>();
    }

    template<class ElementSolution>
    Scalar porosity(const Element& element,
                              const SubControlVolume& scv,
                              const ElementSolution& elemSol) const
    {
        if (pybind11::hasattr(pySpatialParams_, "porosity"))
            return pySpatialParams_.attr("porosity")(element, scv, elemSol).template cast<Scalar>();
        else
            return pySpatialParams_.attr("porosityAtPos")(scv.center()).template cast<Scalar>();
    }

private:
    pybind11::object pySpatialParams_;
};



template <class SP, class... options>
void registerOnePSpatialParams(pybind11::handle scope,
                               pybind11::class_<SP, options...> cls)
{
    using pybind11::operator""_a;
    using GridGeometry = std::decay_t<decltype(std::declval<SP>().gridGeometry())>;

    cls.def(pybind11::init([](std::shared_ptr<GridGeometry> gridGeometry, pybind11::object sp){
        return std::make_shared<SP>(gridGeometry, sp);
    }));

    cls.def("gravity", &SP::gravity);
}

} // namespace Dumux::Python

#endif
