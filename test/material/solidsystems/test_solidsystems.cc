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
 * \ingroup MaterialTests
 * \brief This test makes sure that the programming interface is
 *        observed by all solid systems.
 */

#include <config.h>
#include <iostream>
#include <vector>

#include "checksolidsystem.hh"

// include all solid systems
#include <dumux/material/solidsystems/1csolid.hh>
#include <dumux/material/solidsystems/compositionalsolidphase.hh>

// include all solid components
#include <dumux/material/components/cao.hh>
#include <dumux/material/components/cao2h2.hh>
#include <dumux/material/components/granite.hh>
#include <dumux/material/components/nacl.hh>
#include <dumux/material/components/solid.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;
    using Scalar = double;

    int success = 0;

    ///////////////////////////////////////////////////////////
    // check all solid systems
    // with the solid components CaO, Ca(OH)2, Granite, NaCl
    ///////////////////////////////////////////////////////////

    // check SolidSystems::OneCSolid inert
    {
        using ComponentT = Components::CaO<Scalar>;
        using SolidSystem = SolidSystems::OneCSolid<Scalar, ComponentT, true>;
        success += checkSolidSystem<Scalar, SolidSystem>();
    }
    {
        using ComponentT = Components::CaO2H2<Scalar>;
        using SolidSystem = SolidSystems::OneCSolid<Scalar, ComponentT, true>;
        success += checkSolidSystem<Scalar, SolidSystem>();
    }
    {
        using ComponentT = Components::Granite<Scalar>;
        using SolidSystem = SolidSystems::OneCSolid<Scalar, ComponentT, true>;
        success += checkSolidSystem<Scalar, SolidSystem>();
    }
    {
        using ComponentT = Components::NaCl<Scalar>;
        using SolidSystem = SolidSystems::OneCSolid<Scalar, ComponentT, true>;
        success += checkSolidSystem<Scalar, SolidSystem>();
    }

    // check SolidSystems::OneCSolid non-inert
    {
        using ComponentT = Components::CaO<Scalar>;
        using SolidSystem = SolidSystems::OneCSolid<Scalar, ComponentT, false>;
        success += checkSolidSystem<Scalar, SolidSystem>();
    }
    {
        using ComponentT = Components::CaO2H2<Scalar>;
        using SolidSystem = SolidSystems::OneCSolid<Scalar, ComponentT, false>;
        success += checkSolidSystem<Scalar, SolidSystem>();
    }
    {
        using ComponentT = Components::Granite<Scalar>;
        using SolidSystem = SolidSystems::OneCSolid<Scalar, ComponentT, false>;
        success += checkSolidSystem<Scalar, SolidSystem>();
    }
    {
        using ComponentT = Components::NaCl<Scalar>;
        using SolidSystem = SolidSystems::OneCSolid<Scalar, ComponentT, false>;
        success += checkSolidSystem<Scalar, SolidSystem>();
    }

    // check SolidSystems::CompositionalSolidPhase
    {
        using Component1 = Components::CaO<Scalar>;
        using Component2 = Components::CaO2H2<Scalar>;
        using SolidSystem = SolidSystems::CompositionalSolidPhase<Scalar, Component1, Component2, 2>;
        success += checkSolidSystem<Scalar, SolidSystem>();
    }
    {
        using Component1 = Components::Granite<Scalar>;
        using Component2 = Components::NaCl<Scalar>;
        using SolidSystem = SolidSystems::CompositionalSolidPhase<Scalar, Component1, Component2, 1>;
        success += checkSolidSystem<Scalar, SolidSystem>();
    }
    {
        using Component1 = Components::Granite<Scalar>;
        using Component2 = Components::NaCl<Scalar>;
        using SolidSystem = SolidSystems::CompositionalSolidPhase<Scalar, Component1, Component2, 0>;
        success += checkSolidSystem<Scalar, SolidSystem>();
    }

    return success;

} // end main
