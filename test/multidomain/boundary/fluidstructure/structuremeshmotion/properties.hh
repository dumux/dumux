// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_MULTIDOMAIN_FLUIDSTRUCTURE_STRUCTURE_MESHMOTION_PROPERTIES_HH
#define DUMUX_TEST_MULTIDOMAIN_FLUIDSTRUCTURE_STRUCTURE_MESHMOTION_PROPERTIES_HH

#include <type_traits>

#include <dune/grid/uggrid.hh>
#include <dune/subgrid/subgrid.hh>

#include <dumux/discretization/box.hh>

#include <dumux/geomechanics/hyperelastic/model.hh>
#include <dumux/multidomain/boundary/fluidstructure/meshmotion/model.hh>
#include <dumux/multidomain/boundary/fluidstructure/couplingmanager_structuremesh.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {

struct StructureMeshMotionCommon
{
    using HostGrid = Dune::UGGrid<2>;
    // using HostGrid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
    using Grid = Dune::SubGrid<2, HostGrid>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

struct Structure { using InheritsFrom = std::tuple<StructureMeshMotionCommon, Hyperelastic, BoxModel>; };
// struct MeshMotion { using InheritsFrom = std::tuple<StructureMeshMotionCommon, LinearElasticMeshMotionModel, BoxModel>; };
struct MeshMotion { using InheritsFrom = std::tuple<StructureMeshMotionCommon, BiharmonicMeshMotionModel, BoxModel>; };
} // end namespace TTag

// Set the problem types
template<class TypeTag>
struct Problem<TypeTag, TTag::Structure> { using type = StructureProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Structure>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = DefaultDynamicHyperelasticSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
// struct Problem<TypeTag, TTag::MeshMotion> { using type = LinearElasticMeshMotionProblem<TypeTag>; };
struct Problem<TypeTag, TTag::MeshMotion> { using type = BiharmonicMeshMotionProblem<TypeTag>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::Structure>
{
    using MDTraits = MultiDomainTraits<Properties::TTag::Structure, Properties::TTag::MeshMotion>;
    using type = StructureMeshCouplingManager<MDTraits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::MeshMotion>
{
    using MDTraits = MultiDomainTraits<Properties::TTag::Structure, Properties::TTag::MeshMotion>;
    using type = StructureMeshCouplingManager<MDTraits>;
};

} // end namespace Dumux::Properties

#endif
