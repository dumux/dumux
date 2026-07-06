// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_KIRCHHOFF_LOVE_PLATE_EIGENMODES_PROPERTIES_HH
#define DUMUX_KIRCHHOFF_LOVE_PLATE_EIGENMODES_PROPERTIES_HH

#include <type_traits>

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/box.hh>

#include <dumux/solidmechanics/plate/kirchhoff_love/model.hh>
#include <dumux/solidmechanics/plate/kirchhoff_love/couplingmanager.hh>

#include "../problem.hh"

namespace Dumux::Properties {

namespace TTag {

struct KLPlateEigenmodesCommon
{
    using Grid = Dune::FoamGrid<2, 2>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

struct KLPlateEigenmodesRotation
{ using InheritsFrom = std::tuple<KLPlateEigenmodesCommon, KirchhoffLovePlateRotation, PQ1BubbleModel>; };

struct KLPlateEigenmodesDeformation
{ using InheritsFrom = std::tuple<KLPlateEigenmodesCommon, KirchhoffLovePlateDeformation, BoxModel>; };

} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::KLPlateEigenmodesRotation>
{ using type = KirchhoffLovePlateTestProblemRotation<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::KLPlateEigenmodesDeformation>
{ using type = KirchhoffLovePlateTestProblemDeformation<TypeTag>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::KLPlateEigenmodesRotation>
{
    using MDTraits = MultiDomainTraits<TTag::KLPlateEigenmodesRotation, TTag::KLPlateEigenmodesDeformation>;
    using type = KirchhoffLovePlateCouplingManager<MDTraits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::KLPlateEigenmodesDeformation>
{
    using MDTraits = MultiDomainTraits<TTag::KLPlateEigenmodesRotation, TTag::KLPlateEigenmodesDeformation>;
    using type = KirchhoffLovePlateCouplingManager<MDTraits>;
};

} // end namespace Dumux::Properties

#endif
