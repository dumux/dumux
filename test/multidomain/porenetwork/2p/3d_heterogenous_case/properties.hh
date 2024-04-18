// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief The properties for the two-phase pore network model.
 */
#ifndef DUMUX_PNM2P_PROPERTIES_HH
#define DUMUX_PNM2P_PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porenetwork/2p/model.hh>
#include <dumux/porenetwork/2p/spatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/multishapelocalrules.hh>

#include <dumux/common/properties.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/porenetwork/common/utilities.hh>

#include <dumux/multidomain/porenetwork/constraint/model.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/porenetwork/constraint/couplingmanager.hh>

#include "problem_network.hh"
#include "problem_constraint.hh"
#include "spatialparams_network.hh"

//////////
// Specify the properties
//////////
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DrainageProblem { using InheritsFrom = std::tuple<PNMTwoP>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DrainageProblem> { using type = DrainageProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DrainageProblem>
 {
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Dumux::FluidSystems::H2OAir<Scalar, Dumux::Components::SimpleH2O<Scalar>>;
 };

//! The grid flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::DrainageProblem>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluxVariablesCache = GetPropTypeOr<TypeTag,
        Properties::FluxVariablesCache, FluxVariablesCaching::EmptyCache<Scalar>
    >;
    using Traits = PoreNetwork::PNMTwoPDefaultGridFVCTraits<Problem,
                                                            FluxVariablesCache,
                                                            Dumux::PoreNetwork::TwoPInvasionState<Problem, Dumux::PoreNetwork::StateSwitchMethod::theta>>;
public:
    using type = PoreNetwork::PNMTwoPGridFluxVariablesCache<Problem, FluxVariablesCache, enableCache, Traits>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DrainageProblem>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LocalRules = PoreNetwork::FluidMatrix::MultiShapeTwoPLocalRules<Scalar>;
public:
    using type = PoreNetwork::TwoPDrainageSpatialParams<GridGeometry, Scalar, LocalRules>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DrainageProblem> { using type = Dune::FoamGrid<1, 3>; };

// Create new type tags
namespace TTag {
struct ConstraintProblem { using InheritsFrom = std::tuple<PNMConstraintModel, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ConstraintProblem> { using type = PNMConstraintProblem<TypeTag>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ConstraintProblem> { using type = Dune::FoamGrid<1, 3>; };

// set the coupling manager type
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DrainageProblem>
{
private:
    using Traits = MultiDomainTraits<TTag::DrainageProblem, TTag::ConstraintProblem>;
public:
#if THROATCONSTRAINT
    using type = PNMConstraintCouplingManager< Traits >;
#else
    using type = void;
#endif
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::ConstraintProblem>
{
private:
    using Traits = MultiDomainTraits<TTag::DrainageProblem, TTag::ConstraintProblem>;
public:
#if THROATCONSTRAINT
    using type = PNMConstraintCouplingManager< Traits >;
#else
    using type = void;
#endif
};

} //end namespace Dumux::Properties

#endif
