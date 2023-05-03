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
 * \filegrid
 * \brief Heat conduction problem with multiple solid spheres
 */
#ifndef DUMUX_TEST_PORENETWORK_SOLID_ENERGY_PROPERTIES_HH
#define DUMUX_TEST_PORENETWORK_SOLID_ENERGY_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/porenetwork/solidenergy/model.hh>
#include <dumux/porenetwork/solidenergy/spatialparams.hh>

#include <dumux/material/components/constant.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct PNMSolidModel { using InheritsFrom = std::tuple<PNMSolidEnergy>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMSolidModel>
{ using type = Dumux::SolidProblem<TypeTag>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMSolidModel>
{ using type = Dune::FoamGrid<1, 3>; };

//! The spatial parameters to be employed.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMSolidModel>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::SolidEnergySpatialParams<GridGeometry, Scalar>;
};

// per default the solid system is inert with one constant component
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::PNMSolidModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using InertComponent = Components::Constant<1, Scalar>;
public:
    using type = SolidSystems::InertSolidPhase<Scalar, InertComponent>;
};

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::PNMSolidModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::TruncatedPyramidGrainFouriersLaw<Scalar>; //from grainfourierslaw.hh (specified in solidenergy/model.hh)
};

} // end namespace Dumux::Properties

#endif
