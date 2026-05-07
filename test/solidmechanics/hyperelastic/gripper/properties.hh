// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_GRIPPER_PROPERTIES_HH
#define DUMUX_GRIPPER_PROPERTIES_HH

#include <type_traits>

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/multidomain/traits.hh>

#include "model.hh"
#include "spatialparams.hh"
#include "problem_momentum.hh"
#include "problem_mass.hh"

namespace Dumux::Properties {

// ============================================================
// Momentum type tag: PQ1BubbleModel for displacement
// ============================================================
namespace TTag {
struct GripperMomentum
{
    using InheritsFrom = std::tuple<GripperMomentumModel, PQ1BubbleModel>;

    using Grid = Dune::FoamGrid<2, 2>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;

    template<class TypeTag>
    using Problem = GripperMomentumProblem<TypeTag>;

    template<class TypeTag>
    using SpatialParams = GripperSpatialParams<
        GetPropType<TypeTag, Properties::GridGeometry>,
        GetPropType<TypeTag, Properties::Scalar>
    >;
};
} // end TTag

// ============================================================
// Mass type tag: BoxModel (P1) for pore pressure
// ============================================================
namespace TTag {
struct GripperMass
{
    using InheritsFrom = std::tuple<GripperMassModel, BoxModel>;

    using Grid = Dune::FoamGrid<2, 2>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;

    template<class TypeTag>
    using Problem = GripperMassProblem<TypeTag>;

    template<class TypeTag>
    using SpatialParams = GripperSpatialParams<
        GetPropType<TypeTag, Properties::GridGeometry>,
        GetPropType<TypeTag, Properties::Scalar>
    >;
};
} // end TTag

// ============================================================
// Coupling manager property (set on both type tags)
// ============================================================

} // end namespace Dumux::Properties

namespace Dumux {
//! Convenience alias for the gripper multidomain traits.
using GripperMDTraits = MultiDomainTraits<Properties::TTag::GripperMomentum,
                                          Properties::TTag::GripperMass>;
} // end namespace Dumux

namespace Dumux::Properties {

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::GripperMomentum>
{ using type = GripperCouplingManager<GripperMDTraits>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::GripperMass>
{ using type = GripperCouplingManager<GripperMDTraits>; };

} // end namespace Dumux::Properties

#endif // DUMUX_GRIPPER_PROPERTIES_HH
