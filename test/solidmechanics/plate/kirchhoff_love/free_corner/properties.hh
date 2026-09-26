// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_KIRCHHOFF_LOVE_PLATE_FREE_CORNER_TEST_PROPERTIES_HH
#define DUMUX_KIRCHHOFF_LOVE_PLATE_FREE_CORNER_TEST_PROPERTIES_HH

#include <type_traits>

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/box.hh>

#include <dumux/solidmechanics/plate/kirchhoff_love/model.hh>
#include <dumux/solidmechanics/plate/kirchhoff_love/couplingmanager.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {

struct KLFreeCornerTestCommon
{
    using Grid = Dune::FoamGrid<2, 2>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

struct KLFreeCornerTestRotation
{ using InheritsFrom = std::tuple<KLFreeCornerTestCommon, KirchhoffLovePlateRotation, PQ1BubbleModel>; };

struct KLFreeCornerTestDeformation
{ using InheritsFrom = std::tuple<KLFreeCornerTestCommon, KirchhoffLovePlateDeformation, BoxModel>; };

} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::KLFreeCornerTestRotation>
{ using type = FreeCornerProblemRotation<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::KLFreeCornerTestDeformation>
{ using type = FreeCornerProblemDeformation<TypeTag>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::KLFreeCornerTestRotation>
{
    using MDTraits = MultiDomainTraits<TTag::KLFreeCornerTestRotation, TTag::KLFreeCornerTestDeformation>;
    using type = KirchhoffLovePlateCouplingManager<MDTraits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::KLFreeCornerTestDeformation>
{
    using MDTraits = MultiDomainTraits<TTag::KLFreeCornerTestRotation, TTag::KLFreeCornerTestDeformation>;
    using type = KirchhoffLovePlateCouplingManager<MDTraits>;
};

} // end namespace Dumux::Properties

#endif
