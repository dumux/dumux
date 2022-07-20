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
 *
 * \brief The properties for the two-phase pore network model.
 */
#ifndef DUMUX_PNM2P_PROPERTIES_HH
#define DUMUX_PNM2P_PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porenetwork/2p/model.hh>
#include <dumux/porenetwork/2p/spatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/multishapelocalrules.hh>

#include <dumux/common/properties.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/porenetwork/common/utilities.hh>

#include "problem.hh"
#include "spatialparams.hh"
#include "iofields.hh"

//////////
// Specify the properties
//////////
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
#if ISOTHERMAL
struct DrainageProblem { using InheritsFrom = std::tuple<PNMTwoP>; };
#else
struct DrainageProblem { using InheritsFrom = std::tuple<PNMTwoPNI>; };
#endif
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DrainageProblem> { using type = DrainageProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DrainageProblem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::H2OAir<Scalar, Components::SimpleH2O<Scalar>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DrainageProblem>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LocalRules = PoreNetwork::FluidMatrix::MultiShapeTwoPLocalRules<Scalar>;
public:
    using type = PoreNetwork::TwoPDrainageSpatialParams<GridGeometry, Scalar, LocalRules>;
};

 //!< Set the default formulation to pwsn
template<class TypeTag>
struct Formulation<TypeTag, TTag::DrainageProblem>
{ static constexpr auto value = TwoPFormulation::p0s1; };

#if ISOTHERMAL
template<class TypeTag>
struct IOFields<TypeTag, TTag::DrainageProblem> { using type = PoreNetwork::TwoPIOFieldsGasWet; };
#else
//! Set the vtk output fields specific to the non-isothermal two-phase model
template<class TypeTag>
struct IOFields<TypeTag, TTag::DrainageProblem> { using type = EnergyIOFields<PoreNetwork::TwoPIOFieldsGasWet>; };
#endif

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DrainageProblem> { using type = Dune::FoamGrid<1, DIM>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DrainageProblem> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DrainageProblem> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DrainageProblem> { static constexpr bool value = false; };

} //end namespace Dumux::Properties

#endif
