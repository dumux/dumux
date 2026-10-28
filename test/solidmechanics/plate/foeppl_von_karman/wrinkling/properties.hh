// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_FOEPPL_VON_KARMAN_WRINKLING_TEST_PROPERTIES_HH
#define DUMUX_FOEPPL_VON_KARMAN_WRINKLING_TEST_PROPERTIES_HH

#include <type_traits>

#include <dune/grid/uggrid.hh>

#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/box.hh>

#include <dumux/multidomain/traits.hh>

#include <dumux/solidmechanics/plate/foeppl_von_karman/model.hh>
#include <dumux/solidmechanics/plate/foeppl_von_karman/couplingmanager.hh>
#include <dumux/solidmechanics/plate/foeppl_von_karman/holegauge.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {

struct WrinklingTestCommon
{
    using Grid = Dune::UGGrid<2>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

struct WrinklingTestRotation
{ using InheritsFrom = std::tuple<WrinklingTestCommon, FoepplVonKarmanPlateRotation, PQ1BubbleModel>; };

struct WrinklingTestDeformation
{ using InheritsFrom = std::tuple<WrinklingTestCommon, FoepplVonKarmanPlateDeformation, BoxModel>; };

struct WrinklingTestInPlane
{ using InheritsFrom = std::tuple<WrinklingTestCommon, FoepplVonKarmanPlateInPlane, BoxModel>; };

} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::WrinklingTestRotation>
{ using type = WrinklingTestProblemRotation<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::WrinklingTestDeformation>
{ using type = WrinklingTestProblemDeformation<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::WrinklingTestInPlane>
{ using type = WrinklingTestProblemInPlane<TypeTag>; };

namespace Detail {
using WrinklingTestMDTraits = MultiDomainTraits<TTag::WrinklingTestRotation,
                                                TTag::WrinklingTestDeformation,
                                                TTag::WrinklingTestInPlane>;
} // end namespace Detail

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::WrinklingTestRotation>
{ using type = HoleGaugeCouplingManager<Detail::WrinklingTestMDTraits>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::WrinklingTestDeformation>
{ using type = HoleGaugeCouplingManager<Detail::WrinklingTestMDTraits>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::WrinklingTestInPlane>
{ using type = HoleGaugeCouplingManager<Detail::WrinklingTestMDTraits>; };

} // end namespace Dumux::Properties

#endif
