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
 * \brief The properties for the  one-phase two-component pore network model.
 */
#ifndef DUMUX_PNM1P2C_PROPERTIES_HH
#define DUMUX_PNM1P2C_PROPERTIES_HH

// Pore network model
#include <dumux/porenetworkflow/1pnc/model.hh>

#include <dumux/common/properties.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

// the problem
#include "problem.hh"

//////////
// Specify the properties
//////////
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
#if ISOTHERMAL
struct PNMOnePTwoCProblem { using InheritsFrom = std::tuple<PNMOnePNC>; };
#else
struct PNMOnePTwoCProblem { using InheritsFrom = std::tuple<PNMOnePNCNI>; };
#endif
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMOnePTwoCProblem> { using type = Dumux::PNMOnePTwoCProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMOnePTwoCProblem>
{
    using Policy = FluidSystems::H2ON2DefaultPolicy</*simple*/true>;
    using H2ON2 = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>, Policy>;
    using type = FluidSystems::OnePAdapter<H2ON2>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMOnePTwoCProblem> { using type = Dune::FoamGrid<1, 3>; };

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMOnePTwoCProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = Dumux::PoreNetwork::TransmissibilityAzzamDullien<Scalar>;
public:
    using type = Dumux::PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

} //end namespace Dumux::Properties

#endif
