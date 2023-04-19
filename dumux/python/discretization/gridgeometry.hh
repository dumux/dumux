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

        using Element = typename FVEG::Element;
        cls.def("bind", [](FVEG& self, const Element& element){
            self.bind(element);
        }, "element"_a);
        cls.def("bindElement", [](FVEG& self, const Element& element){
            self.bindElement(element);
        }, "element"_a);

        cls.def_property_readonly("scvs", [](FVEG& self){
            const auto range = scvs(self);
            return pybind11::make_iterator(range.begin(), range.end());
        }, pybind11::keep_alive<0, 1>());
        cls.def_property_readonly("scvfs", [](FVEG& self){
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
        return std::make_shared<GG>(gridView);
    }), "gridView"_a);

    cls.def("update", [](GG& self, const GridView& gridView){
        return self.update(gridView);
    }, "gridView"_a);

    cls.def_property_readonly("numDofs", &GG::numDofs);
    cls.def_property_readonly("numScv", &GG::numScv);
    cls.def_property_readonly("numScvf", &GG::numScvf);
    cls.def_property_readonly("bBoxMax", &GG::bBoxMax);
    cls.def_property_readonly("bBoxMin", &GG::bBoxMin);
    cls.def_property_readonly("gridView", &GG::gridView);

    cls.def_property_readonly_static("discMethod", [](const pybind11::object&){
        return GG::discMethod.name();
    });

    cls.def_property_readonly("localView", [](GG& self){
        return localView(self);
    });

    using LocalView = typename GG::LocalView;
    using Element = typename LocalView::Element;
    cls.def("boundLocalView", [](GG& self, const Element& element){
        auto view = localView(self);
        view.bind(element);
        return view;
    }, "element"_a);

    using SubControlVolume = typename GG::SubControlVolume;
    Impl::registerSubControlVolume<SubControlVolume>(scope);

    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    Impl::registerSubControlVolumeFace<SubControlVolumeFace>(scope);

    // also compile the corresponding local view
    Impl::registerFVElementGeometry<LocalView>(scope);

}

} // namespace Dumux::Python

#endif
