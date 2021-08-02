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
 * \ingroup TwoPTests
 * \brief The properties for the incompressible 2p-boxdfm test.
 */

#ifndef DUMUX_INCOMPRESSIBLE_TWOPBOXDFM_TEST_PROPERTIES_HH
#define DUMUX_INCOMPRESSIBLE_TWOPBOXDFM_TEST_PROPERTIES_HH

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/boxdfm/model.hh>
#include <dumux/porousmediumflow/2p/incompressiblelocalresidual.hh>

#include "problem.hh"
#include "spatialparams.hh"

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2>;
#endif

namespace Dumux::Properties {

// we need to derive first from twop and then from the box-dfm Model
// because the flux variables cache type of TwoP is overwritten in BoxDfmModel
// Create new type tags
namespace TTag {
struct TwoPIncompressibleBoxDfm { using InheritsFrom = std::tuple<BoxDfmModel, TwoP>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPIncompressibleBoxDfm> { using type = GRIDTYPE; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPIncompressibleBoxDfm> { using type = TwoPTestProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPIncompressibleBoxDfm>
{
private:
    using FVG = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = TwoPTestSpatialParams<FVG, Scalar>;
};

// the local residual containing the analytic derivative methods
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::TwoPIncompressibleBoxDfm> { using type = TwoPIncompressibleLocalResidual<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPIncompressibleBoxDfm>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TwoPIncompressibleBoxDfm> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TwoPIncompressibleBoxDfm> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TwoPIncompressibleBoxDfm> { static constexpr bool value = false; };

// Enable the box-interface solver
template<class TypeTag>
struct EnableBoxInterfaceSolver<TypeTag, TTag::TwoPIncompressibleBoxDfm> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
