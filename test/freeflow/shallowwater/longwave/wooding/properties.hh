// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief The properties of the Wooding V-catchment tests
 */
#ifndef DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/freeflow/shallowwater/model.hh>
#include <dumux/freeflow/shallowwater/longwave/model.hh>

#include "spatialparams.hh"
#include "spatialparams_localinertial.hh"
#include "problem.hh"
#include "problem_localinertial.hh"
#include "localinertialflux.hh"

namespace Dumux::Properties::TTag {

struct WoodingLongWave
{
    using InheritsFrom = std::tuple<LongWave, BoxModel>;
    using Grid = Dune::YaspGrid<2, Dune::TensorProductCoordinates<double, 2>>;

    template<class TypeTag>
    using Problem = WoodingProblem<TypeTag>;

    template<class TypeTag>
    using SpatialParams = WoodingSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                               GetPropType<TypeTag, Properties::Scalar>>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

struct WoodingLocalInertial
{
    using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>;
    using Grid = Dune::YaspGrid<2, Dune::TensorProductCoordinates<double, 2>>;

    template<class TypeTag>
    using Problem = WoodingLocalInertialProblem<TypeTag>;

    template<class TypeTag>
    using SpatialParams = WoodingLocalInertialSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                                            GetPropType<TypeTag, Properties::Scalar>,
                                                            GetPropType<TypeTag, Properties::VolumeVariables>>;

    template<class TypeTag>
    using AdvectionType = Wooding::LocalInertialFlux<Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

} // end namespace Dumux::Properties::TTag

#endif
