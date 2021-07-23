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

#ifndef DUMUX_PYTHON_DISCRETIZATION_GRIDGEOMETRY_HH
#define DUMUX_PYTHON_DISCRETIZATION_GRIDGEOMETRY_HH

#include <memory>
#include <dune/common/classname.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/common/typeregistry.hh>

namespace Dumux::Python {

namespace Impl {

template<class SCV>
void registerSubControlVolume(pybind11::handle scope)
{
    using namespace Dune::Python;

    auto [cls, addedToRegistry] = insertClass<SCV>(
        scope, "SubControlVolume",
        GenerateTypeName(Dune::className<SCV>()),
        IncludeFiles{"dumux/python/discretization/gridgeometry.hh"}
    );

    if (addedToRegistry)
    {
        using pybind11::operator""_a;

        cls.def_property_readonly("center", &SCV::center);
        cls.def_property_readonly("volume", &SCV::volume);
        cls.def_property_readonly("dofIndex", &SCV::dofIndex);
        cls.def_property_readonly("localDofIndex", &SCV::localDofIndex);
        cls.def_property_readonly("dofPosition", &SCV::dofPosition);
        cls.def_property_readonly("elementIndex", &SCV::elementIndex);

    }
}

template<class SCVF>
void registerSubControlVolumeFace(pybind11::handle scope)
{
    using namespace Dune::Python;

    auto [cls, addedToRegistry] = insertClass<SCVF>(
        scope, "SubControlVolumeFace",
        GenerateTypeName(Dune::className<SCVF>()),
        IncludeFiles{"dumux/python/discretization/gridgeometry.hh"}
    );

    if (addedToRegistry)
    {
        using pybind11::operator""_a;

        cls.def_property_readonly("center", &SCVF::center);
        cls.def_property_readonly("area", &SCVF::area);
        cls.def_property_readonly("ipGlobal", &SCVF::ipGlobal);
        cls.def_property_readonly("boundary", &SCVF::boundary);
        cls.def_property_readonly("unitOuterNormal", &SCVF::unitOuterNormal);
        cls.def_property_readonly("insideScvIdx", &SCVF::insideScvIdx);
        cls.def_property_readonly("outsideScvIdx", [](SCVF& self){ return self.outsideScvIdx(); });
        cls.def_property_readonly("index", &SCVF::index);
    }
}

template<class FVEG>
void registerFVElementGeometry(pybind11::handle scope)
{
    using namespace Dune::Python;

    auto [cls, addedToRegistry] = insertClass<FVEG>(
        scope, "FVElementGeometry",
        GenerateTypeName(Dune::className<FVEG>()),
        IncludeFiles{"dumux/python/discretization/gridgeometry.hh"}
    );

    if (addedToRegistry)
    {
        using pybind11::operator""_a;

        cls.def_property_readonly("numScvf", &FVEG::numScv);
        cls.def_property_readonly("numScv", &FVEG::numScv);
        cls.def_property_readonly("hasBoundaryScvf", &FVEG::hasBoundaryScvf);
        cls.def("bind", &FVEG::bind, "element"_a);
        cls.def("scvs", [](FVEG& self){
            const auto range = scvs(self);
            return pybind11::make_iterator(range.begin(), range.end());
        }, pybind11::keep_alive<0, 1>());
        cls.def("scvfs", [](FVEG& self){
            const auto range = scvfs(self);
            return pybind11::make_iterator(range.begin(), range.end());
        }, pybind11::keep_alive<0, 1>());
    }
}

} // end namespace Impl

// see python/dumux/discretization/__init__.py for how this is used for JIT compilation
template <class GG, class... Options>
void registerGridGeometry(pybind11::handle scope, pybind11::class_<GG, Options...> cls)
{
    using pybind11::operator""_a;

    using GridView = typename GG::GridView;

    cls.def(pybind11::init([](const GridView& gridView){
        auto gg = std::make_shared<GG>(gridView);
        // remove this once the interface is changed on the C++ side (see #1056)
        gg->update();
        return gg;
    }), "gridView"_a);

    // update this once the interface is changed on the C++ side (see #1056)
    cls.def("update", [](GG& self, const GridView& gridView){
        return self.update();
    }, "gridView"_a);

    cls.def_property_readonly("numDofs", &GG::numDofs);
    cls.def_property_readonly("numScv", &GG::numScv);
    cls.def_property_readonly("numScvf", &GG::numScvf);
    cls.def_property_readonly("bBoxMax", &GG::bBoxMax);
    cls.def_property_readonly("bBoxMin", &GG::bBoxMin);
    cls.def_property_readonly("gridView", &GG::gridView);

    cls.def_property_readonly_static("discMethod", [](const pybind11::object&){
        return toString(GG::discMethod);
    });

    using SubControlVolume = typename GG::SubControlVolume;
    Impl::registerSubControlVolume<SubControlVolume>(scope);

    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    Impl::registerSubControlVolumeFace<SubControlVolumeFace>(scope);

    // also compile the corresponding local view
    using LocalView = typename GG::LocalView;
    Impl::registerFVElementGeometry<LocalView>(scope);
    // and make it accessible
    cls.def("localView", [](GG& self){ return localView(self); });
}

} // namespace Dumux::Python

#endif
