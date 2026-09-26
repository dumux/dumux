// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_KIRCHHOFF_LOVE_PLATE_NAVIER_SQUARE_TEST_PROPERTIES_HH
#define DUMUX_KIRCHHOFF_LOVE_PLATE_NAVIER_SQUARE_TEST_PROPERTIES_HH

#include <type_traits>

#include <dune/grid/uggrid.hh>

#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/box.hh>

#include <dumux/solidmechanics/plate/kirchhoff_love/model.hh>
#include <dumux/solidmechanics/plate/kirchhoff_love/couplingmanager.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {

struct KLNavierSquareTestCommon
{
    using Grid = Dune::UGGrid<2>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

struct KLNavierSquareTestRotation
{ using InheritsFrom = std::tuple<KLNavierSquareTestCommon, KirchhoffLovePlateRotation, PQ1BubbleModel>; };

struct KLNavierSquareTestDeformation
{ using InheritsFrom = std::tuple<KLNavierSquareTestCommon, KirchhoffLovePlateDeformation, BoxModel>; };

} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::KLNavierSquareTestRotation>
{ using type = NavierSquareProblemRotation<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::KLNavierSquareTestDeformation>
{ using type = NavierSquareProblemDeformation<TypeTag>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::KLNavierSquareTestRotation>
{
    using MDTraits = MultiDomainTraits<TTag::KLNavierSquareTestRotation, TTag::KLNavierSquareTestDeformation>;
    using type = KirchhoffLovePlateCouplingManager<MDTraits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::KLNavierSquareTestDeformation>
{
    using MDTraits = MultiDomainTraits<TTag::KLNavierSquareTestRotation, TTag::KLNavierSquareTestDeformation>;
    using type = KirchhoffLovePlateCouplingManager<MDTraits>;
};

} // end namespace Dumux::Properties

#endif
