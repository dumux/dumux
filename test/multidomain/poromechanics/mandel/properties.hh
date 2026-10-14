// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief The properties of Mandel's problem
 */
#ifndef DUMUX_TEST_MANDEL_PROPERTIES_HH
#define DUMUX_TEST_MANDEL_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/poromechanics/poroelastic/model.hh>
#include <dumux/poromechanics/poroelastic/couplingmanager.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/traits.hh>

#include "spatialparams_1p.hh"
#include "spatialparams_poroelastic.hh"
#include "problem_1p.hh"
#include "problem_poroelastic.hh"

namespace Dumux::Properties {

namespace TTag {
struct OnePSub { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
struct PoroElasticSub { using InheritsFrom = std::tuple<PoroElastic, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::OnePSub> { using type = Dune::YaspGrid<2>; };
template<class TypeTag>
struct Grid<TypeTag, TTag::PoroElasticSub> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::OnePSub> { using type = OnePSubProblem<TypeTag>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::PoroElasticSub> { using type = PoroElasticSubProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePSub>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PoroElasticSub>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePSub>
{
    using type = OnePSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                   GetPropType<TypeTag, Properties::Scalar>,
                                   GetPropType<TypeTag, Properties::CouplingManager>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PoroElasticSub>
{
    using type = PoroElasticSpatialParams<GetPropType<TypeTag, Properties::Scalar>,
                                          GetPropType<TypeTag, Properties::GridGeometry>,
                                          GetPropType<TypeTag, Properties::CouplingManager>>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePSub>
{ using type = PoroMechanicsCouplingManager<MultiDomainTraits<TTag::OnePSub, TTag::PoroElasticSub>>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PoroElasticSub>
{ using type = PoroMechanicsCouplingManager<MultiDomainTraits<TTag::OnePSub, TTag::PoroElasticSub>>; };

} // end namespace Dumux::Properties

#endif
