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
 * \ingroup ThreePWaterOilTests
 * \brief The properties of the non-isothermal steam-assisted gravity
 *        drainage (SAGD) problem.
 */

#ifndef DUMUX_SAGDPROPERTIES_HH
#define DUMUX_SAGDPROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/3pwateroil/model.hh>

#include <dumux/material/fluidsystems/h2oheavyoil.hh>
#include <dumux/material/solidsystems/1csolid.hh>
#include <dumux/material/components/constant.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Sagd { using InheritsFrom = std::tuple<ThreePWaterOilNI>; };
struct ThreePWaterOilSagdBox { using InheritsFrom = std::tuple<Sagd, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Sagd> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Sagd> { using type = Dumux::SagdProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Sagd>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = SagdSpatialParams<GridGeometry, Scalar>;
};

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Sagd>
{ using type = Dumux::FluidSystems::H2OHeavyOil<GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag>
struct OnlyGasPhaseCanDisappear<TypeTag, TTag::Sagd> { static constexpr bool value = true; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::Sagd> { static constexpr bool value = true; };

// Set the solid system
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::Sagd>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using InertComponent = Components::Constant<1, Scalar>;
    using type = SolidSystems::InertSolidPhase<Scalar, InertComponent>;
};

} // end namespace Dumux::Properties

#endif
