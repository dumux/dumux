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
 * \brief Properties for for the bulk problem in the 1pnc facet coupling test.
 */
#ifndef DUMUX_TEST_FACETCOUPLING_ONEPNC_BULK_PROPERTIES_HH
#define DUMUX_TEST_FACETCOUPLING_ONEPNC_BULK_PROPERTIES_HH

#ifndef DIMWORLD
#define DIMWORLD 2
#endif

#include <dune/alugrid/grid.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include "problem_bulk.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// create the type tag nodes
namespace TTag {
struct BaseBulk {};
struct OnePNCBulk { using InheritsFrom = std::tuple<OnePNC, BaseBulk>; };
struct OnePNCNIBulk { using InheritsFrom = std::tuple<OnePNCNI, BaseBulk>; };

struct OnePNCBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePNCBulk>; };
struct OnePNCNIBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePNCNIBulk>; };

struct OnePNCBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePNCBulk>; };
struct OnePNCNIBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePNCNIBulk>; };

struct OnePNCBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, OnePNCBulk>; };
struct OnePNCNIBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, OnePNCNIBulk>; };
} // end namespace TTag

// Set the grid type (DIMWORLD is defined in CMakeLists.txt)
template<class TypeTag>
struct Grid<TypeTag, TTag::BaseBulk> { using type = Dune::ALUGrid<2, DIMWORLD, Dune::cube, Dune::nonconforming>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::BaseBulk> { using type = OnePNCBulkProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BaseBulk>
{
    using type = OnePSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BaseBulk>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2ON2 = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*simplified=*/true>>;
public:
    using type = FluidSystems::OnePAdapter<H2ON2, H2ON2::liquidPhaseIdx>;
};

} // end namespace Dumux::Properties

#endif
