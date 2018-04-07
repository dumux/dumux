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
 * \brief This test makes sure that the programming interface is
 *        observed by all fluid systems
 */
#include <config.h>

#include "checkfluidsystem.hh"

// include all fluid systems in dumux-stable
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/fluidsystems/brineair.hh>
#include <dumux/material/fluidsystems/brineco2.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidsystems/h2oairmesitylene.hh>
#include <dumux/material/fluidsystems/h2oairxylene.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/h2on2kinetic.hh>
#include <dumux/material/fluidsystems/h2on2o2.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/spe5.hh>

// include all fluid states
#include <dumux/material/fluidstates/2p2c.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/fluidstates/isothermalimmiscible.hh>
#include <dumux/material/fluidstates/nonequilibrium.hh>
#include <dumux/material/fluidstates/nonequilibriumenergy.hh>
#include <dumux/material/fluidstates/nonequilibriummass.hh>
#include <dumux/material/fluidstates/pressureoverlay.hh>
#include <dumux/material/fluidstates/pseudo1p2c.hh>
#include <dumux/material/fluidstates/saturationoverlay.hh>
#include <dumux/material/fluidstates/temperatureoverlay.hh>

int main()
{
    using namespace Dumux;

    using Scalar = double;
    using H2O = Components::H2O<Scalar>;
    using N2 = Components::N2<Scalar>;

    using Liquid = FluidSystems::OnePLiquid<Scalar, H2O>;
    using Gas = FluidSystems::OnePGas<Scalar, N2>;

    int success = 0;
    std::vector<std::string> collectedExceptions;

    /////////////////////////
    // check all fluid states
    {
        using FluidSystem = FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false>;
        using BaseFluidState = CompositionalFluidState<Scalar, FluidSystem>;
        BaseFluidState baseFs;

        // TwoPTwoCFluidState TODO: doesn't fulfill interface!
        // TwoPTwoCFluidState<Scalar, FluidSystem> fluidStateTwoPTwoC;
        // success += checkFluidState<Scalar>(fluidStateTwoPTwoC);

        // CompositionalFluidState
        CompositionalFluidState<Scalar, FluidSystem> compositionalFluidState;
        success += checkFluidState<Scalar>(compositionalFluidState);

        // ImmiscibleFluidState
        ImmiscibleFluidState<Scalar, FluidSystem> immiscibleFluidState;
        success += checkFluidState<Scalar>(immiscibleFluidState);

        // IsothermalImmiscibleFluidState
        IsothermalImmiscibleFluidState<Scalar, FluidSystem> isothermalImmiscibleFluidState;
        success += checkFluidState<Scalar>(isothermalImmiscibleFluidState);

        // NonEquilibriumFluidState
        NonEquilibriumFluidState<Scalar, FluidSystem> nonEquilibriumFluidState;
        success += checkFluidState<Scalar>(nonEquilibriumFluidState);

        // NonEquilibriumEnergyFluidState TODO: fails the test
        // NonEquilibriumEnergyFluidState<Scalar, FluidSystem> nonEquilibriumEnergyFluidState;
        // success += checkFluidState<Scalar>(nonEquilibriumEnergyFluidState);

        // NonEquilibriumMassFluidState
        NonEquilibriumMassFluidState<Scalar, FluidSystem> nonEquilibriumMassFluidState;
        success += checkFluidState<Scalar>(nonEquilibriumMassFluidState);

        // PressureOverlayFluidState
        PressureOverlayFluidState<Scalar, BaseFluidState> pressureOverlayFluidState(baseFs);
        success += checkFluidState<Scalar>(pressureOverlayFluidState);

        // SaturationOverlayFluidState
        SaturationOverlayFluidState<Scalar, BaseFluidState> saturationOverlayFluidState(baseFs);
        success += checkFluidState<Scalar>(saturationOverlayFluidState);

        // TemperatureOverlayFluidState
        TemperatureOverlayFluidState<Scalar, BaseFluidState> temperatureOverlayFluidState(baseFs);
        success += checkFluidState<Scalar>(temperatureOverlayFluidState);
    }


    //////////////////////////
    // check all fluid systems

    // 2p-immiscible
    {   using FluidSystem = FluidSystems::TwoPImmiscible<Scalar, Liquid, Liquid>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = FluidSystems::TwoPImmiscible<Scalar, Liquid, Gas>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = FluidSystems::TwoPImmiscible<Scalar, Gas, Liquid>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // base

    // Brine -- Air
    {   using H2OType = Components::SimpleH2O<Scalar>;
        const bool enableComplexRelations=false;
        using FluidSystem = FluidSystems::BrineAir<Scalar, H2OType, enableComplexRelations>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType = Components::SimpleH2O<Scalar>;
        const bool enableComplexRelations=true;
        using FluidSystem = FluidSystems::BrineAir<Scalar, H2OType, enableComplexRelations>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // Brine -- CO2

    // H2O -- Air
    {   using H2OType = Components::SimpleH2O<Scalar>;
        const bool enableComplexRelations=false;
        using FluidSystem = FluidSystems::H2OAir<Scalar, H2OType, enableComplexRelations>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType = Components::SimpleH2O<Scalar>;
        const bool enableComplexRelations=true;
        using FluidSystem = FluidSystems::H2OAir<Scalar, H2OType, enableComplexRelations>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType =  Components::H2O<Scalar>;
        const bool enableComplexRelations=false;
        using FluidSystem = FluidSystems::H2OAir<Scalar, H2OType, enableComplexRelations>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType =  Components::H2O<Scalar>;
        const bool enableComplexRelations=true;
        using FluidSystem = FluidSystems::H2OAir<Scalar, H2OType, enableComplexRelations>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType = Components::TabulatedComponent<Components::H2O<Scalar>>;
        const bool enableComplexRelations=false;
        using FluidSystem = FluidSystems::H2OAir<Scalar, H2OType, enableComplexRelations>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType = Components::TabulatedComponent<Components::H2O<Scalar>>;
        const bool enableComplexRelations=true;
        using FluidSystem = FluidSystems::H2OAir<Scalar, H2OType, enableComplexRelations>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // gas phase
    {   using FluidSystem = FluidSystems::OnePGas<Scalar, H2O>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Mesitylene
    {   using FluidSystem = FluidSystems::H2OAirMesitylene<Scalar>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Xylene
    {   using FluidSystem = FluidSystems::H2OAirXylene<Scalar>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2
    {   using FluidSystem = FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/true>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2 -- kinetic
    {   using FluidSystem = FluidSystems::H2ON2Kinetic<Scalar, /*enableComplexRelations=*/false>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = FluidSystems::H2ON2Kinetic<Scalar, /*enableComplexRelations=*/true>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2 -- o2
    {   using FluidSystem = FluidSystems::H2ON2O2<Scalar, /*enableComplexRelations=*/false>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = FluidSystems::H2ON2O2<Scalar, /*enableComplexRelations=*/true>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // liquid phase
    {   using FluidSystem = FluidSystems::OnePLiquid<Scalar, H2O>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // spe5
    {   using FluidSystem = FluidSystems::Spe5<Scalar>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    return success;
}
