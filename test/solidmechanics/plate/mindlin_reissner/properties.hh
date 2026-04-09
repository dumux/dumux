// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_MINDLIN_REISSNER_PLATE_TEST_PROPERTIES_HH
#define DUMUX_MINDLIN_REISSNER_PLATE_TEST_PROPERTIES_HH

#include <type_traits>

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/box.hh>

#include <dumux/solidmechanics/plate/mindlin_reissner/model.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {

struct MRPlateTestCommon
{
    using Grid = Dune::FoamGrid<2, 2>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

struct MRPlateTestRotation
{ using InheritsFrom = std::tuple<MRPlateTestCommon, MindlinReissnerPlateRotation, PQ1BubbleModel>; };

struct MRPlateTestDeformation
{ using InheritsFrom = std::tuple<MRPlateTestCommon, MindlinReissnerPlateDeformation, BoxModel>; };

} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::MRPlateTestRotation>
{ using type = MindlinReissnerPlateTestProblemRotation<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::MRPlateTestDeformation>
{ using type = MindlinReissnerPlateTestProblemDeformation<TypeTag>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::MRPlateTestRotation>
{
    using MDTraits = MultiDomainTraits<TTag::MRPlateTestRotation, TTag::MRPlateTestDeformation>;
    using type = MindlinReissnerPlateCouplingManager<MDTraits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::MRPlateTestDeformation>
{
    using MDTraits = MultiDomainTraits<TTag::MRPlateTestRotation, TTag::MRPlateTestDeformation>;
    using type = MindlinReissnerPlateCouplingManager<MDTraits>;
};

} // end namespace Dumux::Properties

#endif
