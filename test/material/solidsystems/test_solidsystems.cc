// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief This test makes sure that the programming interface is
 *        observed by all solid systems.
 */

#include <config.h>
#include <cmath>
#include <iostream>
#include <utility>
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

#include <dumux/material/solidstates/compositionalsolidstate.hh>

namespace Dumux {

//! Returns 1 if the heat capacity times the density of the solid phase is not
//! the volume average of the heat capacities per volume of its components
template<class Scalar, class SolidSystem>
int checkCompositionalHeatCapacity(Scalar volFrac0, Scalar volFrac1)
{
    using Component0 = typename SolidSystem::ComponentOne;
    using Component1 = typename SolidSystem::ComponentTwo;

    CompositionalSolidState<Scalar, SolidSystem> solidState;
    const Scalar temperature = 293.15;
    solidState.setTemperature(temperature);
    solidState.setVolumeFraction(SolidSystem::comp0Idx, volFrac0);
    solidState.setVolumeFraction(SolidSystem::comp1Idx, volFrac1);

    const Scalar expected = volFrac0*Component0::solidDensity(temperature)*Component0::solidHeatCapacity(temperature)
                            + volFrac1*Component1::solidDensity(temperature)*Component1::solidHeatCapacity(temperature);
    const Scalar actual = (volFrac0 + volFrac1)*SolidSystem::density(solidState)*SolidSystem::heatCapacity(solidState);

    if (std::abs(actual - expected) < 1e-12*expected)
        return 0;

    std::cout << "heat capacity per volume of " << SolidSystem::componentName(SolidSystem::comp0Idx)
              << " (volume fraction " << volFrac0 << ") and " << SolidSystem::componentName(SolidSystem::comp1Idx)
              << " (volume fraction " << volFrac1 << "): " << actual << " instead of " << expected << std::endl;
    return 1;
}

} // end namespace Dumux

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

    // check the heat capacity of SolidSystems::CompositionalSolidPhase
    for (const auto& [volFrac0, volFrac1] : std::vector<std::pair<Scalar, Scalar>>{{0.6, 0.1}, {0.05, 0.7}, {0.95, 0.005}})
    {
        using CaOSystem = SolidSystems::CompositionalSolidPhase<Scalar, Components::CaO<Scalar>, Components::CaO2H2<Scalar>>;
        success += checkCompositionalHeatCapacity<Scalar, CaOSystem>(volFrac0, volFrac1);
        using GraniteSystem = SolidSystems::CompositionalSolidPhase<Scalar, Components::Granite<Scalar>, Components::NaCl<Scalar>>;
        success += checkCompositionalHeatCapacity<Scalar, GraniteSystem>(volFrac0, volFrac1);
    }

    return success;

} // end main
