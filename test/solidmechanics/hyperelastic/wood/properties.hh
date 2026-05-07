// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_WOOD_TEST_PROPERTIES_HH
#define DUMUX_WOOD_TEST_PROPERTIES_HH

#include <type_traits>

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/box.hh>

#include "model.hh"
#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties::TTag {

struct WoodTest
{
    using InheritsFrom = std::tuple<Wood, BoxModel>;

    using Grid = Dune::YaspGrid<2>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;

    template<class TypeTag>
    using Problem = WoodProblem<TypeTag>;

    template<class TypeTag>
    using SpatialParams = WoodSpatialParams<
        GetPropType<TypeTag, Properties::GridGeometry>,
        GetPropType<TypeTag, Properties::Scalar>
    >;
};

} // end namespace Dumux::Properties

#endif
