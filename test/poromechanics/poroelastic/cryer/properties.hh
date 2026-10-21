// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Properties of the small-deformation Cryer benchmark: one-phase flow (CC-MPFA)
 *        coupled to linear poroelasticity (box).
 */
#ifndef DUMUX_TEST_POROELASTIC_CRYER_PROPERTIES_HH
#define DUMUX_TEST_POROELASTIC_CRYER_PROPERTIES_HH

#include <tuple>

#include <dune/alugrid/grid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/ccmpfa.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/poromechanics/poroelastic/model.hh>

#include "couplingmanager.hh"

#include "problem_1p.hh"
#include "problem_poroelastic.hh"
#include "spatialparams_1p.hh"
#include "spatialparams_poroelastic.hh"

namespace Dumux::Properties {

namespace TTag {
struct CryerOneP { using InheritsFrom = std::tuple<OneP, CCMpfaModel>; };
struct CryerPoroElastic { using InheritsFrom = std::tuple<PoroElastic, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::CryerOneP> { using type = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>; };
template<class TypeTag>
struct Grid<TypeTag, TTag::CryerPoroElastic> { using type = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::CryerOneP> { using type = CryerOnePProblem<TypeTag>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::CryerPoroElastic> { using type = CryerPoroElasticProblem<TypeTag>; };

// permeability and fluid density do not depend on the solution, so the MPFA transmissibilities are
// assembled once and kept in the grid flux variables cache
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::CryerOneP> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::CryerOneP> { static constexpr bool value = true; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::CryerOneP>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::CryerPoroElastic>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CryerOneP>
{
    using type = CryerOnePSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                        GetPropType<TypeTag, Properties::Scalar>,
                                        GetPropType<TypeTag, Properties::CouplingManager>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CryerPoroElastic>
{
    using type = CryerPoroElasticSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                               GetPropType<TypeTag, Properties::Scalar>,
                                               GetPropType<TypeTag, Properties::CouplingManager>>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::CryerOneP>
{ using type = CryerCouplingManager<MultiDomainTraits<TTag::CryerOneP, TTag::CryerPoroElastic>>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::CryerPoroElastic>
{ using type = CryerCouplingManager<MultiDomainTraits<TTag::CryerOneP, TTag::CryerPoroElastic>>; };

} // end namespace Dumux::Properties

#endif
