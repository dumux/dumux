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

#ifndef DUMUX_MATERIAL_COMPONENTS_SUSPENDEDBIOMASS_HH
#define DUMUX_MATERIAL_COMPONENTS_SUSPENDEDBIOMASS_HH

// ## The suspended biomass component (`suspendedbiomass.hh`)
//
// This file contains the __ component class__ which defines the name and molar mass of suspended biomass
//
// [[content]]
//
// ### Include files
// [[codeblock]]
// including the base component
#include <dumux/material/components/base.hh>

#include <dumux/common/parameters.hh>
// [[/codeblock]]

// ### The suspended biomass component
// [[codeblock]]
namespace Dumux::Components {

// In SuspendedBiomass, we define the properties of the component suspended biomass
template <class Scalar>
class SuspendedBiomass
: public Components::Base<Scalar, SuspendedBiomass<Scalar> >
{
public:
    // the name
    static std::string name()
    { return "Suspended_Biomass"; }
// [[/codeblock]]

// ### The suspended biomass component's properties
// [[codeblock]]
    // The molar mass, which is not really defined for suspended biomass. Thus, we read it from params.input or use a default of 1.
    // Based on a cell mass of 2.5e-16, the molar mass of cells would be 1.5e8 kg/mol, but such high molar masses would lead to numerical problems.
    static Scalar molarMass()
    {
        Scalar molarMass = getParam<Scalar>("BioCoefficients.SuspendedBiomassMolarMass", 1);
        return molarMass;
    }
};


} // end namespace Dumux::Components
// [[/codeblock]]
// [[/content]]
#endif
