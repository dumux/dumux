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

#ifndef DUMUX_PYTHON_POROUSMEDIUMFLOW_FVSPATIALPARAMS_1P_HH
#define DUMUX_PYTHON_POROUSMEDIUMFLOW_FVSPATIALPARAMS_1P_HH

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <dumux/python/common/fvspatialparams.hh>

namespace Dumux::Python {

template<class GridGeometry_, class PT>
class FVSpatialParamsOneP
: public Dumux::Python::FVSpatialParams<GridGeometry_>
{
    using ThisType = Dumux::Python::FVSpatialParamsOneP<GridGeometry_, PT>;
    using ParentType = Dumux::Python::FVSpatialParams<GridGeometry_>;
public:
    using GridGeometry = GridGeometry_;
    using Scalar = typename GridGeometry::GridView::ctype;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using PermeabilityType = PT;

    FVSpatialParamsOneP(std::shared_ptr<const GridGeometry> gridGeometry,
                        pybind11::object pySpatialParams)
    : ParentType(gridGeometry, pySpatialParams)
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

    template<class SolidSystem, class ElementSolution>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const

    {
        if (pybind11::hasattr(pySpatialParams_, "inertVolumeFraction"))
            return pySpatialParams_.attr("inertVolumeFraction")(element, scv, elemSol, compIdx).template cast<Scalar>();
        else if (pybind11::hasattr(pySpatialParams_, "inertVolumeFractionAtPos"))
            return pySpatialParams_.attr("inertVolumeFractionAtPos")(scv.center(), compIdx).template cast<Scalar>();
        else
            return 1.0 - this->porosity(element, scv, elemSol);
    }

    static constexpr bool evaluatePermeabilityAtScvfIP()
    { return false; }

private:
    pybind11::object pySpatialParams_;
};

template <class SpatialParams, class... options>
void registerFVSpatialParamsOneP(pybind11::handle scope, pybind11::class_<SpatialParams, options...> cls)
{
    using pybind11::operator""_a;
    using GridGeometry = typename SpatialParams::GridGeometry;

    cls.def(pybind11::init([](std::shared_ptr<const GridGeometry> gridGeometry, pybind11::object p){
        return std::make_shared<SpatialParams>(gridGeometry, p);
    }));

    cls.def("permeability", &SpatialParams::template permeability<decltype(std::ignore)>);
    cls.def("porosity", &SpatialParams::template porosity<decltype(std::ignore)>);
    cls.def("inertVolumeFraction", &SpatialParams::template inertVolumeFraction<decltype(std::ignore), decltype(std::ignore)>);
}

} // namespace Dumux::Python

#endif
