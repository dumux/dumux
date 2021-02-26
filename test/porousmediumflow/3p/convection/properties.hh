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
 * \ingroup ThreePTests
 * \brief The properties of the test for the ThreePModel in combination with the NI model for a convection problem.
 */

#ifndef DUMUX_3PNI_CONVECTION_PROPERTIES_HH
#define DUMUX_3PNI_CONVECTION_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/3p/model.hh>
#include <dumux/material/fluidsystems/h2oairmesitylene.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidmatrixinteractions/3p/thermalconductivitysomerton3p.hh>

#include "../conduction/spatialparams.hh"  // reuse the conduction spatialParams
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct ThreePNIConvection { using InheritsFrom = std::tuple<ThreePNI>; };
struct ThreePNIConvectionBox { using InheritsFrom = std::tuple<ThreePNIConvection, BoxModel>; };
struct ThreePNIConvectionCCTpfa { using InheritsFrom = std::tuple<ThreePNIConvection, CCTpfaModel>; };
struct ThreePNIConvectionCCMpfa { using InheritsFrom = std::tuple<ThreePNIConvection, CCMpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ThreePNIConvection> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ThreePNIConvection> { using type = ThreePNIConvectionProblem<TypeTag>; };


// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ThreePNIConvection>
{ using type = FluidSystems::H2OAirMesitylene<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::ThreePNIConvection>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ThreePNISpatialParams<GridGeometry, Scalar>;
};

} // end namespace Dumux::Properties

#endif
