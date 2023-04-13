// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_MATERIAL_COMPONENTS_BIOFILM_HH
#define DUMUX_MATERIAL_COMPONENTS_BIOFILM_HH

// ## The biofilm component (`biofilm.hh`)
//
// This file contains the __solid component class__ which defines the name, molar mass and density of biofilm
//
// [[content]]
//
// ### Include files
// [[codeblock]]
// including the base and the generic solid component
#include <dumux/material/components/base.hh>
#include <dumux/material/components/solid.hh>

#include <dumux/common/parameters.hh>
// [[/codeblock]]

// ### The biofilm component
// [[codeblock]]
namespace Dumux::Components {

// In Biofilm, we define the properties of the solid component biofilm
template <class Scalar>
class Biofilm
: public Components::Base<Scalar, Biofilm<Scalar> >
, public Components::Solid<Scalar, Biofilm<Scalar> >
{
public:
    // the name
    static std::string name()
    { return "Biofilm"; }
// [[/codeblock]]

// ### The biofilm component's properties
// [[codeblock]]
    // The molar mass, which is not really defined for biofilm. Thus, we read it from params.input or use a default of 1.
    // Based on a cell mass of 2.5e-16, the molar mass of cells would be 1.5e8 kg/mol, but biofilms are more than just cells and such high molar masses would lead to numerical problems.
    static Scalar molarMass()
    {
        static Scalar molarMass = getParam<Scalar>("BioCoefficients.BiofilmMolarMass", 1);
        return molarMass;
    }

    // The density, or rather the dry density (dry biomass/wet volume), most of biofilm is water.
    // It is typically highly variable for different biofilms, thus we read it from params.input or use the default value of 10 kg/m^3
    static Scalar solidDensity(Scalar temperature)
    {
        static Scalar rho = getParam<Scalar>("BioCoefficients.BiofilmDensity", 10);
        return rho;
    }
};

} // end namespace Dumux::Components
// [[/codeblock]]
// [[/content]]
#endif
