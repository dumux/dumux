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
 * \ingroup NavierStokesTests
 * \brief Test for the staggered grid Navier-Stokes model with analytical solution.
 */
#ifndef DUMUX_KOVASZNAY_TEST_PROPERTIES_MOMENTUM_HH
#define DUMUX_KOVASZNAY_TEST_PROPERTIES_MOMENTUM_HH

#ifndef ENABLECACHING
#define ENABLECACHING true
#endif

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2>
#endif

#ifndef NAVIER_STOKES_MODEL
#define NAVIER_STOKES_MODEL NavierStokesMomentum
#endif

#ifndef DISCRETIZATION_MODEL
#define DISCRETIZATION_MODEL FaceCenteredStaggeredModel
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/flux/fluxvariablescaching.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/discretization/fcstaggered.hh>

#include <dumux/freeflow/navierstokes/momentum/diamond/model.hh>
#include <dumux/discretization/fcdiamond.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct KovasznayTestMomentum { using InheritsFrom = std::tuple<NAVIER_STOKES_MODEL, DISCRETIZATION_MODEL>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::KovasznayTestMomentum>
{
    using type = Dumux::KovasznayTestProblem<TypeTag>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::KovasznayTestMomentum>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::KovasznayTestMomentum>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = std::conditional_t<
        GridGeometry::discMethod == DiscretizationMethods::fcdiamond,
        FaceCenteredDiamondFluxVariablesCache<Scalar, GridGeometry>,
        FluxVariablesCaching::EmptyCache<Scalar>
    >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::KovasznayTestMomentum> { using type = GRIDTYPE; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::KovasznayTestMomentum> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::KovasznayTestMomentum> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::KovasznayTestMomentum> { static constexpr bool value = ENABLECACHING; };

} // end namespace Dumux::Properties

#endif
