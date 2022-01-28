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
 * \ingroup PoromechanicsTests
 * \brief Definition of the spatial parameters for the single-phase flow
 *        sub-problem in the coupled poro-mechanical el1p problem.
 */
#ifndef DUMUX_GEOMECHANICS_ELONEP_PROPERTIES_HH
#define DUMUX_GEOMECHANICS_ELONEP_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dumux/io/grid/gridmanager_alu.hh>
#endif
#if HAVE_DUNE_UGGRID
#include <dumux/io/grid/gridmanager_ug.hh>
#endif

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/geomechanics/poroelastic/model.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/geomechanics/poroelastic/couplingmanager.hh>

#include "spatialparams_1p.hh"
#include "spatialparams_poroelastic.hh"

#include "problem_1p.hh"
#include "problem_poroelastic.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct OnePSub { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
} // end namespace TTag

// The fluid phase consists of one constant component
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePSub>
{
    using type = Dumux::FluidSystems::OnePLiquid< GetPropType<TypeTag, Properties::Scalar>,
                                                  Dumux::Components::Constant<0, GetPropType<TypeTag, Properties::Scalar>> >;
};

// Set the grid type
template<class TypeTag>
#ifndef GRIDTYPE
struct Grid<TypeTag, TTag::OnePSub> { using type = Dune::YaspGrid<2>; };
#else
struct Grid<TypeTag, TTag::OnePSub> { using type = GRIDTYPE; };
#endif
// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePSub> { using type = OnePSubProblem<TypeTag> ; };
// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePSub>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using type = OnePSpatialParams<GridGeometry, Scalar, CouplingManager>;
};

// Create new type tags
namespace TTag {
struct PoroElasticSub { using InheritsFrom = std::tuple<PoroElastic, BoxModel>; };
} // end namespace TTag
// Set the grid type
template<class TypeTag>
#ifndef GRIDTYPE
struct Grid<TypeTag, TTag::PoroElasticSub> { using type = Dune::YaspGrid<2>; };
#else
struct Grid<TypeTag, TTag::PoroElasticSub> { using type = GRIDTYPE; };
#endif
// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PoroElasticSub> { using type = Dumux::PoroElasticSubProblem<TypeTag>; };
// The fluid phase consists of one constant component
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PoroElasticSub>
{
    using type = Dumux::FluidSystems::OnePLiquid< GetPropType<TypeTag, Properties::Scalar>,
                                                  Dumux::Components::Constant<0, GetPropType<TypeTag, Properties::Scalar>> >;
};
// The spatial parameters property
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PoroElasticSub>
{
    using type = PoroElasticSpatialParams< GetPropType<TypeTag, Properties::Scalar>,
                                           GetPropType<TypeTag, Properties::GridGeometry> >;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePSub>
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<TTag::OnePSub, TTag::PoroElasticSub>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PoroElasticSub>
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<TTag::OnePSub, TTag::PoroElasticSub>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
};

} // end namespace Dumux::Properties

#endif
