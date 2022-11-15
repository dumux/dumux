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

#ifndef DUMUX_MICP_SUB_DOMAIN_PROPERTIES_HH
#define DUMUX_MICP_SUB_DOMAIN_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2pncmin/model.hh>


// We include the necessary material files
#include <examples/biomineralization/material/fluidsystems/biominsimplechemistry.hh>
#include <examples/biomineralization/material/solidsystems/biominsolids.hh>
#include <examples/biomineralization/material/co2tables.hh>

// We include the problem and spatial parameters headers used for this simulation.
#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// We create new type tag for our simulation which inherits from the 2pncmin model and the box discretization
namespace TTag {
struct MICPColumnSimpleChemistry { using InheritsFrom = std::tuple<TwoPNCMin, BoxModel>; };
} // end namespace TTag
// [[/codeblock]]

// ### Property specializations
//
// In the following piece of code, mandatory `properties` for which no meaningful
// default can be set, are specialized for our type tag `MICPColumnSimpleChemistry`.

// [[codeblock]]
// We set the grid to a 1D Yasp Grid
template<class TypeTag>
struct Grid<TypeTag, TTag::MICPColumnSimpleChemistry>
{ using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };

// We set the problem  used for our simulation, defining boundary and initial conditions (see below)
template<class TypeTag>
struct Problem<TypeTag, TTag::MICPColumnSimpleChemistry>
{ using type = MICPColumnProblemSimpleChemistry<TypeTag>; };

// We set the fluidSystem  used for our simulation
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MICPColumnSimpleChemistry>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CO2Tables = BiomineralizationCO2Tables::CO2Tables;
    using H2OTabulated = Components::TabulatedComponent<Components::H2O<Scalar>>;
    using type = Dumux::FluidSystems::BioMinSimpleChemistryFluid<Scalar, CO2Tables, H2OTabulated>;
};

// We set the solidSystem  used for our simulation
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::MICPColumnSimpleChemistry>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = SolidSystems::BioMinSolidPhase<Scalar>;
};

// We define the spatial parameters for our simulation. The values are specified in the corresponding spatialparameters header file, which is included above.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::MICPColumnSimpleChemistry>
{
    using type = ICPSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                  GetPropType<TypeTag, Properties::Scalar>>;
};

// We set the two-phase primary variable formulation used for our simulation
template<class TypeTag>
struct Formulation<TypeTag, TTag::MICPColumnSimpleChemistry>
{ static constexpr auto value = TwoPFormulation::p0s1; };

}// We leave the namespace Dumux::Properties.
// [[/codeblock]]
// [[/content]]
#endif
