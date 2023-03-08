// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Base class for all models which use the one-phase,
 *        fully implicit model.
 *        Adaption of the fully implicit scheme to the one-phase pore network model.
 *
 * NEW: Adaption to solidenergy based on PorousMediumFlow-SolidEnergy
 */

#ifndef DUMUX_PNM_SOLID_ENERGY_MODEL_HH
#define DUMUX_PNM_SOLID_ENERGY_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/porenetwork/properties.hh>

#include <dumux/porousmediumflow/solidenergy/model.hh>
#include <dumux/flux/porenetwork/grainfourierslaw.hh>

#include "spatialparams.hh"
#include "volumevariables.hh"
#include "fluxvariablescache.hh"
#include "iofields.hh"

/*!
 * \ingroup Pore network model
 * \brief A one-phase, isothermal flow model using the fully implicit scheme.
 *
 * one-phase, isothermal flow model, which uses a standard Darcy approach as the
 * equation for the conservation of momentum:
 * \f[
 v = - \frac{\textbf K}{\mu}
 \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 * \f]
 *
 * and solves the mass continuity equation:
 * \f[
 \phi \frac{\partial \varrho}{\partial t} + \text{div} \left\lbrace
 - \varrho \frac{\textbf K}{\mu} \left( \textbf{grad}\, p -\varrho {\textbf g} \right) \right\rbrace = q,
 * \f]
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 * The model supports compressible as well as incompressible fluids.
 */

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Dumux::Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit single-phase problems
// Create new type tags
namespace TTag {
struct PNMSolidEnergy{ using InheritsFrom = std::tuple<PoreNetworkModel, SolidEnergy>; };
}

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PNMSolidEnergy>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    using Traits = SolidEnergyVolumeVariablesTraits<PV, SSY, SST, MT>;
public:
    using type = PoreNetwork::SolidEnergyVolumeVariables<Traits>;
};

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::PNMSolidEnergy>
{ using type = PoreNetwork::TruncatedPyramidGrainFouriersLaw<GetPropType<TypeTag, Properties::Scalar>>; };

//! The spatial parameters to be employed.
//! Use PNMOnePSpatialParams by default.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMSolidEnergy>
{
private:
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::SolidEnergySpatialParams<FVGridGeometry, Scalar>;
};

//! The flux variables cache
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PNMSolidEnergy>
{ using type = PoreNetwork::SolidEnergyFluxVariablesCache<GetPropType<TypeTag, Properties::Scalar>>; };

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMSolidEnergy>
{ using type = PoreNetwork::SolidEnergyIOFields; };

} // namespace Dumux::Properies

#endif
