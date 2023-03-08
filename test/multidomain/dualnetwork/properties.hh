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
 * \brief Heat problem with multiple solid spheres
 */
#ifndef DUMUX_TEST_MULTIDOMAIN_DUALNETWORK_PROPERTIES_HH
#define DUMUX_TEST_MULTIDOMAIN_DUALNETWORK_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager_sub.hh>

#include <dumux/porenetwork/solidenergy/model.hh>

#include "problem_solid.hh"
#include "spatialparams.hh"
#include "fourierslaw.hh"
#include "constant.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct PNMSolidModel { using InheritsFrom = std::tuple<PNMSolidEnergy>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMSolidModel>
{ using type = Dumux::SolidSubProblem<TypeTag>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMSolidModel>
{ using type = Dune::SubGrid<1, Dune::FoamGrid<1, 3>>; };

//! The spatial parameters to be employed.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMSolidModel>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = SolidSpatialParams<GridGeometry, Scalar>;
};

// per default the solid system is inert with one constant component
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::PNMSolidModel>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using InertComponent = Components::NewConstant<1, Scalar>;
    using type = SolidSystems::InertSolidPhase<Scalar, InertComponent>;
};

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::PNMSolidModel>
{ using type = PoreNetwork::FlexibleFouriersLaw<false>; };


} // end namespace Dumux::Properties

#include <dumux/porenetwork/1p/model.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include "problem_void.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct PNMVoidModel { using InheritsFrom = std::tuple<PNMOnePNI>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMVoidModel>
{ using type = Dumux::VoidSubProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMVoidModel>
{
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::NewConstant<2, Scalar>> ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMVoidModel>
{ using type = Dune::SubGrid<1, Dune::FoamGrid<1, 3>>; };

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMVoidModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = PoreNetwork::TransmissibilityPatzekSilin<Scalar, false/*considerPoreBodyResistance*/>;
public:
    using type = PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::PNMVoidModel>
{ using type = PoreNetwork::FlexibleFouriersLaw<true>; };

// use the incompressible local residual (provides analytic jacobian)
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PNMVoidModel>
{ using type = OnePIncompressibleLocalResidual<TypeTag>; };

//! The spatial parameters to be employed.
//! Use PNMOnePSpatialParams by default.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMVoidModel>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Dumux::Properties

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/dualnetwork/couplingmanager.hh>

namespace Dumux::Properties {

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PNMSolidModel>
{
    using Traits = MultiDomainTraits<TypeTag, Properties::TTag::PNMVoidModel>;
    using type = PNMHeatTransferCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PNMVoidModel>
{
    using Traits = MultiDomainTraits<Properties::TTag::PNMSolidModel, TypeTag>;
    using type = PNMHeatTransferCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
