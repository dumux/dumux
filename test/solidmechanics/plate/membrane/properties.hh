// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_MEMBRANE_PLATE_TEST_PROPERTIES_HH
#define DUMUX_MEMBRANE_PLATE_TEST_PROPERTIES_HH

#include <type_traits>

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/discretization/box.hh>

#include <dumux/solidmechanics/plate/membrane/model.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct MembranePlateTest
{
    using InheritsFrom = std::tuple<MembranePlate, BoxModel>;
    using Grid = Dune::FoamGrid<2, 2>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::MembranePlateTest>
{ using type = MembranePlateTestProblem<TypeTag>; };

} // end namespace Dumux::Properties

#endif
