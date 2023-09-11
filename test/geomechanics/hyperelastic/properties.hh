// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTICITY_TEST_PROPERTIES_HH
#define DUMUX_HYPERELASTICITY_TEST_PROPERTIES_HH

#include <type_traits>

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/box.hh>
#include <dumux/geomechanics/hyperelastic/model.hh>

#include "problem.hh"

namespace Dumux::Properties::TTag {

struct HyperelasticityTest
{
    using InheritsFrom = std::tuple<Hyperelastic, BoxModel>;

    using Grid = Dune::YaspGrid<3>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;

    template<class TypeTag>
    using Problem = HyperelasticityProblem<TypeTag>;
};

} // end namespace Dumux::Properties

#endif
