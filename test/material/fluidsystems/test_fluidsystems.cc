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
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidsystems/h2oairmesitylene.hh>
#include <dumux/material/fluidsystems/h2oairxylene.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/h2on2kinetic.hh>
#include <dumux/material/fluidsystems/h2on2o2.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/fluidsystems/purewatersimple.hh>
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
    using Scalar = double;
    using H2O = Dumux::H2O<Scalar>;
    using N2 = Dumux::N2<Scalar>;

    using Liquid = Dumux::FluidSystems::LiquidPhase<Scalar, H2O>;
    using Gas = Dumux::FluidSystems::GasPhase<Scalar, N2>;

    int success = 0;
    std::vector<std::string> collectedExceptions;

    /////////////////////////
    // check all fluid states
    {
        using FluidSystem = Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false>;
        using BaseFluidState = Dumux::CompositionalFluidState<Scalar, FluidSystem>;
        BaseFluidState baseFs;

        // TwoPTwoCFluidState TODO: doesn't fulfill interface!
        // Dumux::TwoPTwoCFluidState<Scalar, FluidSystem> fluidStateTwoPTwoC;
        // success += Dumux::checkFluidState<Scalar>(fluidStateTwoPTwoC);

        // CompositionalFluidState
        Dumux::CompositionalFluidState<Scalar, FluidSystem> compositionalFluidState;
        success += Dumux::checkFluidState<Scalar>(compositionalFluidState);

        // ImmiscibleFluidState
        Dumux::ImmiscibleFluidState<Scalar, FluidSystem> immiscibleFluidState;
        success += Dumux::checkFluidState<Scalar>(immiscibleFluidState);

        // IsothermalImmiscibleFluidState
        Dumux::IsothermalImmiscibleFluidState<Scalar, FluidSystem> isothermalImmiscibleFluidState;
        success += Dumux::checkFluidState<Scalar>(isothermalImmiscibleFluidState);

        // NonEquilibriumFluidState
        Dumux::NonEquilibriumFluidState<Scalar, FluidSystem> nonEquilibriumFluidState;
        success += Dumux::checkFluidState<Scalar>(nonEquilibriumFluidState);

        // NonEquilibriumEnergyFluidState TODO: fails the test
        // Dumux::NonEquilibriumEnergyFluidState<Scalar, FluidSystem> nonEquilibriumEnergyFluidState;
        // success += Dumux::checkFluidState<Scalar>(nonEquilibriumEnergyFluidState);

        // NonEquilibriumMassFluidState
        Dumux::NonEquilibriumMassFluidState<Scalar, FluidSystem> nonEquilibriumMassFluidState;
        success += Dumux::checkFluidState<Scalar>(nonEquilibriumMassFluidState);

        // PressureOverlayFluidState
        Dumux::PressureOverlayFluidState<Scalar, BaseFluidState> pressureOverlayFluidState(baseFs);
        success += Dumux::checkFluidState<Scalar>(pressureOverlayFluidState);

        // SaturationOverlayFluidState
        Dumux::SaturationOverlayFluidState<Scalar, BaseFluidState> saturationOverlayFluidState(baseFs);
        success += Dumux::checkFluidState<Scalar>(saturationOverlayFluidState);

        // TemperatureOverlayFluidState
        Dumux::TemperatureOverlayFluidState<Scalar, BaseFluidState> temperatureOverlayFluidState(baseFs);
        success += Dumux::checkFluidState<Scalar>(temperatureOverlayFluidState);
    }


    //////////////////////////
    // check all fluid systems

    // 2p-immiscible
    {   using FluidSystem = Dumux::FluidSystems::TwoPImmiscible<Scalar, Liquid, Liquid>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = Dumux::FluidSystems::TwoPImmiscible<Scalar, Liquid, Gas>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = Dumux::FluidSystems::TwoPImmiscible<Scalar, Gas, Liquid>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // base

    // Brine -- Air
    {   using H2OType =  Dumux::SimpleH2O<Scalar>;
        const bool enableComplexRelations=false;
        using FluidSystem = Dumux::FluidSystems::BrineAir<Scalar, H2OType, enableComplexRelations>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType =  Dumux::SimpleH2O<Scalar>;
        const bool enableComplexRelations=true;
        using FluidSystem = Dumux::FluidSystems::BrineAir<Scalar, H2OType, enableComplexRelations>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // Brine -- CO2

    // H2O -- Air
    {   using H2OType =  Dumux::SimpleH2O<Scalar>;
        const bool enableComplexRelations=false;
        using FluidSystem = Dumux::FluidSystems::H2OAir<Scalar, H2OType, enableComplexRelations>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType =  Dumux::SimpleH2O<Scalar>;
        const bool enableComplexRelations=true;
        using FluidSystem = Dumux::FluidSystems::H2OAir<Scalar, H2OType, enableComplexRelations>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType =  Dumux::H2O<Scalar>;
        const bool enableComplexRelations=false;
        using FluidSystem = Dumux::FluidSystems::H2OAir<Scalar, H2OType, enableComplexRelations>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType =  Dumux::H2O<Scalar>;
        const bool enableComplexRelations=true;
        using FluidSystem = Dumux::FluidSystems::H2OAir<Scalar, H2OType, enableComplexRelations>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType = Dumux::TabulatedComponent<Scalar, Dumux::H2O<Scalar>>;
        const bool enableComplexRelations=false;
        using FluidSystem = Dumux::FluidSystems::H2OAir<Scalar, H2OType, enableComplexRelations>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType = Dumux::TabulatedComponent<Scalar, Dumux::H2O<Scalar>>;
        const bool enableComplexRelations=true;
        using FluidSystem = Dumux::FluidSystems::H2OAir<Scalar, H2OType, enableComplexRelations>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // gas phase
    {   using FluidSystem = Dumux::FluidSystems::GasPhase<Scalar, H2O>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Mesitylene
    {   using FluidSystem = Dumux::FluidSystems::H2OAirMesitylene<Scalar>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Xylene
    {   using FluidSystem = Dumux::FluidSystems::H2OAirXylene<Scalar>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2
    {   using FluidSystem = Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/true>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2 -- kinetic
    {   using FluidSystem = Dumux::FluidSystems::H2ON2Kinetic<Scalar, /*enableComplexRelations=*/false>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = Dumux::FluidSystems::H2ON2Kinetic<Scalar, /*enableComplexRelations=*/true>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2 -- o2
    {   using FluidSystem = Dumux::FluidSystems::H2ON2O2<Scalar, /*enableComplexRelations=*/false>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = Dumux::FluidSystems::H2ON2O2<Scalar, /*enableComplexRelations=*/true>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // liquid phase
    {   using FluidSystem = Dumux::FluidSystems::LiquidPhase<Scalar, H2O>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // pure water simple
    {   using FluidSystem = Dumux::FluidSystems::PureWaterSimpleFluidSystem<Scalar, /*enableComplexRelations=*/false>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = Dumux::FluidSystems::PureWaterSimpleFluidSystem<Scalar, /*enableComplexRelations=*/true>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // spe5
    {   using FluidSystem = Dumux::FluidSystems::Spe5<Scalar>;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    return success;
}
