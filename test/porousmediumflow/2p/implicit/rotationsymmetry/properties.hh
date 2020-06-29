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

#ifndef DUMUX_TEST_TWOP_ROTATIONALSYMMETRY_PROPERTIES_HH
#define DUMUX_TEST_TWOP_ROTATIONALSYMMETRY_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/common/properties.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/porousmediumflow/2p/model.hh>

#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct TwoPRotationalSymmetryDome
{ using InheritsFrom = std::tuple<TwoP, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPRotationalSymmetryDome>
{ using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPRotationalSymmetryDome>
{ using type = TwoPRotationalSymmetryProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPRotationalSymmetryDome>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::TabulatedComponent<Components::H2O<Scalar>, false>>;
    using NonwettingPhase = FluidSystems::OnePGas<Scalar, Components::TabulatedComponent<Components::Air<Scalar>, false>>;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPRotationalSymmetryDome>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TwoPRotationalSymmetrySpatialParams<GridGeometry, Scalar>;
};

// Set the rotational symmetric grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::TwoPRotationalSymmetryDome>
{
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    struct Traits : public BoxDefaultGridGeometryTraits<GridView> { using Extrusion = RotationalExtrusion<0>; };
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, Traits>;
};

// Caching options
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TwoPRotationalSymmetryDome> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TwoPRotationalSymmetryDome> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TwoPRotationalSymmetryDome> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
