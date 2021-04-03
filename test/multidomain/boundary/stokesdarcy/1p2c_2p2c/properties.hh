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
 * \ingroup BoundaryTests
 * \brief A simple Darcy test problem (cell-centered finite volume method).
 */

#ifndef DUMUX_STOKESDARCY_PROPERTIES_HH
#define DUMUX_STOKESDARCY_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/porousmediumflow/2p2c/model.hh>
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>

#include <dumux/multidomain/boundary/stokesdarcy/couplingmanager.hh>

#include "problem_darcy.hh"
#include "problem_stokes.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
#if !NONISOTHERMAL
struct DarcyTwoPTwoC { using InheritsFrom = std::tuple<TwoPTwoC, CCTpfaModel>; };
#else
struct DarcyTwoPTwoC { using InheritsFrom = std::tuple<TwoPTwoCNI, CCTpfaModel>; };
#endif
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyTwoPTwoC> { using type = Dumux::DarcySubProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyTwoPTwoC> { using type = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>; };

//! Set the default formulation to pw-Sn: This can be over written in the problem.
template<class TypeTag>
struct Formulation<TypeTag, TTag::DarcyTwoPTwoC>
{ static constexpr auto value = TwoPFormulation::p1s0; };

//// The gas component balance (air) is replaced by the total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::DarcyTwoPTwoC> { static constexpr int value = 3; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyTwoPTwoC> { using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::DarcyTwoPTwoC> { static constexpr bool value = true; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyTwoPTwoC>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TwoPTwoCSpatialParams<GridGeometry, Scalar>;
};

// Create new type tags
namespace TTag {
#if !NONISOTHERMAL
struct StokesOnePTwoC { using InheritsFrom = std::tuple<NavierStokesNC, StaggeredFreeFlowModel>; };
#else
struct StokesOnePTwoC { using InheritsFrom = std::tuple<NavierStokesNCNI, StaggeredFreeFlowModel>; };
#endif
} // end namespace TTag


// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::StokesOnePTwoC> { using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// The fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::StokesOnePTwoC>
{
  using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
  static constexpr auto phaseIdx = H2OAir::gasPhaseIdx; // simulate the water phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::StokesOnePTwoC> { static constexpr int value = 3; };

// Use formulation based on mass fractions
template<class TypeTag>
struct UseMoles<TypeTag, TTag::StokesOnePTwoC> { static constexpr bool value = true; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesOnePTwoC> { using type = Dumux::StokesSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::StokesOnePTwoC> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::StokesOnePTwoC> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::StokesOnePTwoC> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::StokesOnePTwoC>
{
    using Traits = StaggeredMultiDomainTraits<TypeTag, TypeTag, Properties::TTag::DarcyTwoPTwoC>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DarcyTwoPTwoC>
{
    using Traits = StaggeredMultiDomainTraits<Properties::TTag::StokesOnePTwoC, Properties::TTag::StokesOnePTwoC, TypeTag>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
