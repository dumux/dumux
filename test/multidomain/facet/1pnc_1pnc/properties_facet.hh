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
 * \brief Properties for the lower-dimensional domain in the 1pnc facet coupling test.
 */
#ifndef DUMUX_TEST_FACETCOUPLING_ONEPNC_FACET_PROPERTIES_HH
#define DUMUX_TEST_FACETCOUPLING_ONEPNC_FACET_PROPERTIES_HH

#ifndef DIMWORLD
#define DIMWORLD 2
#endif

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include "problem_facet.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// create the type tag nodes
namespace TTag {
struct BaseFacet {};
struct OnePNCFacet { using InheritsFrom = std::tuple<OnePNC, BaseFacet>; }; // Isothermal case
struct OnePNCNIFacet { using InheritsFrom = std::tuple<OnePNCNI, BaseFacet>; }; // Non-Isothermal case

struct OnePNCFacetTpfa { using InheritsFrom = std::tuple<OnePNCFacet, CCTpfaModel>; };
struct OnePNCNIFacetTpfa { using InheritsFrom = std::tuple<OnePNCNIFacet, CCTpfaModel>; };

struct OnePNCFacetBox { using InheritsFrom = std::tuple<OnePNCFacet, BoxModel>; };
struct OnePNCNIFacetBox { using InheritsFrom = std::tuple<OnePNCNIFacet, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::BaseFacet> { using type = Dune::FoamGrid<1, DIMWORLD>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::BaseFacet> { using type = OnePNCLowDimProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BaseFacet>
{
    using type = OnePSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BaseFacet>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2ON2 = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*simplified=*/true>>;
public:
    using type = FluidSystems::OnePAdapter<H2ON2, H2ON2::liquidPhaseIdx>;
};

} // end namespace Dumux::Properties

#endif
