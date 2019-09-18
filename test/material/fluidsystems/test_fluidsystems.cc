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
 *        observed by all fluid systems.
 */

#include <config.h>

#include "checkfluidsystem.hh"
#include <dumux/common/parameters.hh>

// include all fluid systems in dumux-stable
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2p1c.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidsystems/3pimmiscible.hh>
#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/fluidsystems/brine.hh>
#include <dumux/material/fluidsystems/brineair.hh>
#include <dumux/material/fluidsystems/brineco2.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidsystems/h2oairmesitylene.hh>
#include <dumux/material/fluidsystems/h2oairxylene.hh>
#include <dumux/material/fluidsystems/h2oheavyoil.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/h2on2kinetic.hh>
#include <dumux/material/fluidsystems/h2on2o2.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>
#include <dumux/material/fluidsystems/spe5.hh>

// include all fluid states
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/fluidstates/isothermalimmiscible.hh>
#include <dumux/material/fluidstates/nonequilibrium.hh>
#include <dumux/material/fluidstates/nonequilibriummass.hh>
#include <dumux/material/fluidstates/pressureoverlay.hh>
#include <dumux/material/fluidstates/pseudo1p2c.hh>
#include <dumux/material/fluidstates/saturationoverlay.hh>
#include <dumux/material/fluidstates/temperatureoverlay.hh>

// for co2, include the tables of the co2 test
#include <test/porousmediumflow/co2/implicit/co2tables.hh>

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
        using FluidSystem = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
        using BaseFluidState = CompositionalFluidState<Scalar, FluidSystem>;
        BaseFluidState baseFs;

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

        // NonEquilibriumMassFluidState
        NonEquilibriumMassFluidState<Scalar, FluidSystem> nonEquilibriumMassFluidState;
        success += checkFluidState<Scalar>(nonEquilibriumMassFluidState);

        // PressureOverlayFluidState
        PressureOverlayFluidState<BaseFluidState> pressureOverlayFluidState(baseFs);
        success += checkFluidState<Scalar>(pressureOverlayFluidState);

        // SaturationOverlayFluidState
        SaturationOverlayFluidState<BaseFluidState> saturationOverlayFluidState(baseFs);
        success += checkFluidState<Scalar>(saturationOverlayFluidState);

        // TemperatureOverlayFluidState
        TemperatureOverlayFluidState<BaseFluidState> temperatureOverlayFluidState(baseFs);
        success += checkFluidState<Scalar>(temperatureOverlayFluidState);
    }


    //////////////////////////
    // check all fluid systems

    // 2p1c
    {   using FluidSystem = FluidSystems::TwoPOneC<Scalar, H2O >;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // 2p-immiscible
    {   using FluidSystem = FluidSystems::TwoPImmiscible<Scalar, Liquid, Liquid>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = FluidSystems::TwoPImmiscible<Scalar, Liquid, Gas>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = FluidSystems::TwoPImmiscible<Scalar, Gas, Liquid>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // 3p-immiscible
    {   using FluidSystem = FluidSystems::ThreePImmiscible<Scalar, Liquid, Liquid, Gas>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // base

    // Brine
    {   using H2OType = Components::SimpleH2O<Scalar>;
        using FluidSystem = FluidSystems::Brine<Scalar, H2OType>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // Brine -- Air
    {   using H2OType = Components::SimpleH2O<Scalar>;
        using FluidSystem = FluidSystems::BrineAir<Scalar, H2OType, FluidSystems::BrineAirDefaultPolicy</*fastButSimplifiedRelations=*/true>>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType = Components::SimpleH2O<Scalar>;
        using FluidSystem = FluidSystems::BrineAir<Scalar, H2OType>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // Brine -- CO2
    // BrineCO2 does not fulfill the restrictToPhase-assertion where we assume that for all
    // functions depending on a phase index only fluid properties of this phase are used
    // that is why checkFluidSystem() needs to be called with "false" here.
    // Also see the checkFluidSystem documentation.
    {   using H2OType = Components::SimpleH2O<Scalar>;
        using FluidSystem = FluidSystems::BrineCO2< Scalar, HeterogeneousCO2Tables::CO2Tables,
                                                    H2OType, FluidSystems::BrineCO2DefaultPolicy</*useConstantSalinity=*/true> >;
        Parameters::init([](auto& params){ params["Brine.Salinity"] = "0.3"; });
        success += checkFluidSystem<Scalar, FluidSystem>( false ); }
    {   using H2OType = Components::SimpleH2O<Scalar>;
        using FluidSystem = FluidSystems::BrineCO2< Scalar, HeterogeneousCO2Tables::CO2Tables,
                                                    H2OType, FluidSystems::BrineCO2DefaultPolicy</*useConstantSalinity=*/false> >;
        success += checkFluidSystem<Scalar, FluidSystem>( false ); }
    {   using H2OType = Components::H2O<Scalar>;
        using FluidSystem = FluidSystems::BrineCO2< Scalar, HeterogeneousCO2Tables::CO2Tables,
                                                    H2OType, FluidSystems::BrineCO2DefaultPolicy</*useConstantSalinity=*/true> >;
        Parameters::init([](auto& params){ params["Brine.Salinity"] = "0.3"; });
        success += checkFluidSystem<Scalar, FluidSystem>( false ); }
    {   using H2OType = Components::H2O<Scalar>;
        using FluidSystem = FluidSystems::BrineCO2< Scalar, HeterogeneousCO2Tables::CO2Tables,
                                                    H2OType, FluidSystems::BrineCO2DefaultPolicy</*useConstantSalinity=*/false> >;
        success += checkFluidSystem<Scalar, FluidSystem>( false ); }
    {   using H2OType = Components::TabulatedComponent<Components::H2O<Scalar>>;
        using FluidSystem = FluidSystems::BrineCO2< Scalar, HeterogeneousCO2Tables::CO2Tables,
                                                    H2OType, FluidSystems::BrineCO2DefaultPolicy</*useConstantSalinity=*/true> >;
        Parameters::init([](auto& params){ params["Brine.Salinity"] = "0.3"; });
        success += checkFluidSystem<Scalar, FluidSystem>( false ); }
    {   using H2OType = Components::TabulatedComponent<Components::H2O<Scalar>>;
        using FluidSystem = FluidSystems::BrineCO2< Scalar, HeterogeneousCO2Tables::CO2Tables,
                                                    H2OType, FluidSystems::BrineCO2DefaultPolicy</*useConstantSalinity=*/false> >;
        success += checkFluidSystem<Scalar, FluidSystem>( false ); }
    {   using H2OType = Components::TabulatedComponent<Components::H2O<Scalar>>;
        using FluidSystem = FluidSystems::BrineCO2< Scalar, HeterogeneousCO2Tables::CO2Tables,
                                                    H2OType, FluidSystems::BrineCO2DefaultPolicy</*useConstantSalinity=*/true,
                                                    /*fastButSimplifiedRelations*/true> >;
        Parameters::init([](auto& params){ params["Brine.Salinity"] = "0.3"; });
        success += checkFluidSystem<Scalar, FluidSystem>( false ); }
    {   using H2OType = Components::TabulatedComponent<Components::H2O<Scalar>>;
        using FluidSystem = FluidSystems::BrineCO2< Scalar, HeterogeneousCO2Tables::CO2Tables,
                                                    H2OType, FluidSystems::BrineCO2DefaultPolicy</*useConstantSalinity=*/false,
                                                    /*fastButSimplifiedRelations*/true> >;
        success += checkFluidSystem<Scalar, FluidSystem>( false ); }

    // H2O -- Air
    {   using H2OType = Components::SimpleH2O<Scalar>;
        using FluidSystem = FluidSystems::H2OAir<Scalar, H2OType, FluidSystems::H2OAirDefaultPolicy</*fastButSimplifiedRelations=*/true>>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType = Components::SimpleH2O<Scalar>;
        using FluidSystem = FluidSystems::H2OAir<Scalar, H2OType>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType =  Components::H2O<Scalar>;
        using FluidSystem = FluidSystems::H2OAir<Scalar, H2OType, FluidSystems::H2OAirDefaultPolicy</*fastButSimplifiedRelations=*/true>>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType =  Components::H2O<Scalar>;
        using FluidSystem = FluidSystems::H2OAir<Scalar, H2OType>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType = Components::TabulatedComponent<Components::H2O<Scalar>>;
        using FluidSystem = FluidSystems::H2OAir<Scalar, H2OType, FluidSystems::H2OAirDefaultPolicy</*fastButSimplifiedRelations=*/true>>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType = Components::TabulatedComponent<Components::H2O<Scalar>>;
        using FluidSystem = FluidSystems::H2OAir<Scalar, H2OType>;
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

    // H2O -- Heavyoil
    {   using FluidSystem = FluidSystems::H2OHeavyOil<Scalar>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2
    {   using FluidSystem = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = FluidSystems::H2ON2<Scalar>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2 -- kinetic
    {   using FluidSystem = FluidSystems::H2ON2Kinetic<Scalar, FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = FluidSystems::H2ON2Kinetic<Scalar>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2 -- o2
    {   using FluidSystem = FluidSystems::H2ON2O2<Scalar, FluidSystems::H2ON2O2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = FluidSystems::H2ON2O2<Scalar>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // liquid phase
    {   using FluidSystem = FluidSystems::OnePLiquid<Scalar, H2O>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // liquid phase 2c
    {
        using FluidSystem = FluidSystems::LiquidPhaseTwoC<Scalar, H2O, Components::Constant<1, Scalar>>;
            Parameters::init([](auto& params){ params["Component.MolarMass"] = "1.0";});
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // spe5
    {   using FluidSystem = FluidSystems::Spe5<Scalar>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    // 1p adapter
    {   using FluidSystem = FluidSystems::OnePAdapter<FluidSystems::TwoPImmiscible<Scalar, Gas, Liquid>, 0>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = FluidSystems::OnePAdapter<FluidSystems::TwoPImmiscible<Scalar, Gas, Liquid>, 1>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using H2OType = Components::TabulatedComponent<Components::H2O<Scalar>>;
        using FluidSystem = FluidSystems::OnePAdapter<FluidSystems::H2OAir<Scalar, H2OType>, 0>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }
    {   using FluidSystem = FluidSystems::OnePAdapter<FluidSystems::H2ON2O2<Scalar>, 1>;
        success += checkFluidSystem<Scalar, FluidSystem>(); }

    return success;
}
