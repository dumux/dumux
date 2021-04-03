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
 * \brief The properties for a simple Navier-Stokes test for the staggered grid (Navier-)Stokes
 * model.
 */
#ifndef DUMUX_STOKES_SUB_PROPERTIES_HH
#define DUMUX_STOKES_SUB_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include "problem_stokes.hh"

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct StokesOnePTwoC { using InheritsFrom = std::tuple<NavierStokesNC, StaggeredFreeFlowModel>; };
} // end namespace TTag

// The fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::StokesOnePTwoC>
{
  using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
  static constexpr auto phaseIdx = H2OAir::liquidPhaseIdx; // simulate the water phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::StokesOnePTwoC> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesOnePTwoC> { using type = Dumux::StokesSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::StokesOnePTwoC> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::StokesOnePTwoC> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::StokesOnePTwoC> { static constexpr bool value = true; };

// Use moles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::StokesOnePTwoC> { static constexpr bool value = true; };

// Do not replace one equation with a total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::StokesOnePTwoC> { static constexpr int value = 3; };
} // end namespace Dumux::Properties

#endif
