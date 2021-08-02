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
 * \brief The problem for the bulk domain in the single-phase facet coupling test.
 */

#ifndef DUMUX_TEST_TPFAFACETCOUPLING_PROPERTIES_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_PROPERTIES_HH

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams.hh"
#include "problem_lowdim.hh"
#include "problem_bulk.hh"

namespace Dumux::Properties {

// create the type tag nodes
// Create new type tags
namespace TTag {
struct OnePBulk { using InheritsFrom = std::tuple<OneP>; };
struct OnePBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePBulk>; };
struct OnePBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, OnePBulk>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePBulk> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };
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

// Create new type tags
namespace TTag {
struct OnePLowDim { using InheritsFrom = std::tuple<OneP>; };
struct OnePLowDimTpfa { using InheritsFrom = std::tuple<OnePLowDim, CCTpfaModel>; };
struct OnePLowDimBox { using InheritsFrom = std::tuple<OnePLowDim, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePLowDim> { using type = Dune::FoamGrid<1, 2>; };
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

// obtain/define some types to be used below in the property definitions and in main
template< class BulkTypeTag, class LowDimTypeTag >
class TestTraits
{
    using BulkFVGridGeometry = Dumux::GetPropType<BulkTypeTag, Dumux::Properties::GridGeometry>;
    using LowDimFVGridGeometry = Dumux::GetPropType<LowDimTypeTag, Dumux::Properties::GridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingMapper<BulkFVGridGeometry, LowDimFVGridGeometry>;
    using CouplingManager = Dumux::FacetCouplingManager<MDTraits, CouplingMapper>;
};

// set cm property for the box test
using BoxTraits = TestTraits<Properties::TTag::OnePBulkBox, Properties::TTag::OnePLowDimBox>;
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkBox> { using type = typename BoxTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePLowDimBox> { using type = typename BoxTraits::CouplingManager; };

// set cm property for the tpfa test
using TpfaTraits = TestTraits<Properties::TTag::OnePBulkTpfa, Properties::TTag::OnePLowDimTpfa>;
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePLowDimTpfa> { using type = typename TpfaTraits::CouplingManager; };

} // end namespace Dumux::Properties

#endif
