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
 * \ingroup RichardsTests
 * \brief The properties of the water infiltration problem with a
 *        low-permeability lens embedded into a high-permeability domain which
 *        uses the Richards box model.
 */
#ifndef DUMUX_RICHARDS_LENSPROPERTIES_HH
#define DUMUX_RICHARDS_LENSPROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "spatialparams.hh"
#include "problem.hh"

// Specify the properties for the lens problem
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RichardsLens { using InheritsFrom = std::tuple<Richards>; };
struct RichardsLensBox { using InheritsFrom = std::tuple<RichardsLens, BoxModel>; };
struct RichardsLensCC { using InheritsFrom = std::tuple<RichardsLens, CCTpfaModel>; };
struct RichardsLensCCMpfa { using InheritsFrom = std::tuple<RichardsLens, CCMpfaModel>; };
} // end namespace TTag

#ifndef GRIDTYPE
// Use 2d YaspGrid
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsLens> { using type = Dune::YaspGrid<2>; };
#else
// Use GRIDTYPE from CMakeLists.txt
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsLens> { using type = GRIDTYPE; };
#endif

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsLens> { using type = RichardsLensProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsLens>
{
    using type = RichardsLensSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                           GetPropType<TypeTag, Properties::Scalar>>;
};

} // end namespace Dumux::Properties

#endif
