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
 * \brief The properties for the bulk domain in the single-phase facet coupling test.
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEP_BULK_PROPERTIES_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEP_BULK_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/properties.hh>

#include "spatialparams.hh"
#include "problem_bulk.hh"

// default for the bulk grid type
#ifndef BULKGRIDTYPE
#define BULKGRIDTYPE Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>
#endif

namespace Dumux::Properties {

// create the type tag nodes
namespace TTag {
struct OnePBulk { using InheritsFrom = std::tuple<OneP>; };
struct OnePBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePBulk>; };
struct OnePBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePBulk>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePBulk> { using type = BULKGRIDTYPE; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePBulk> { using type = OnePBulkProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePBulk>
{
    using type = OnePSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePBulk>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
};

} // end namespace Dumux::Properties

#endif
