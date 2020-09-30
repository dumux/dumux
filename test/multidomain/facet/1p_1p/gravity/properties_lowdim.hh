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
 * \ingroup FacetTests
 * \brief The properties for the lower-dimensional domain in the single-phase facet coupling test.
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEP_LOWDIM_PROPERTIES_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEP_LOWDIM_PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include "spatialparams.hh"
#include "problem_lowdim.hh"

// default for the grid type
#ifndef LOWDIMGRIDTYPE
#define LOWDIMGRIDTYPE Dune::FoamGrid<1, 2>
#endif

namespace Dumux::Properties {

// create the type tag nodes
namespace TTag {
struct OnePLowDim { using InheritsFrom = std::tuple<OneP>; };
struct OnePLowDimTpfa { using InheritsFrom = std::tuple<OnePLowDim, CCTpfaModel>; };

// we need an additional type tag for the test using mpfa in the bulk domain
struct OnePLowDimMpfa { using InheritsFrom = std::tuple<OnePLowDim, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePLowDim> { using type = LOWDIMGRIDTYPE; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePLowDim> { using type = OnePLowDimProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePLowDim>
{
    using type = OnePSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePLowDim>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
};

} // end namespace Dumux::Properties

#endif
