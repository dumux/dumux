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

// include the property system just to make sure that all fluid system
// type tag adapter behave nicely together
#include <dumux/common/propertysystem.hh>

#include "checkfluidsystem.hh"

// include all fluid systems in dumux-stable
#include <dumux/material/fluidsystems/1p.hh>
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
    typedef double Scalar;
    typedef Dumux::H2O<Scalar> H2O;
    typedef Dumux::N2<Scalar> N2;

    typedef Dumux::FluidSystems::LiquidPhase<Scalar, H2O> Liquid;
    typedef Dumux::FluidSystems::GasPhase<Scalar, N2> Gas;

    int success = 0;
    std::string collectedExceptions;

    /////////////////////////
    // check all fluid states
    {
        typedef Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> BaseFluidState;
        BaseFluidState baseFs;

        // TwoPTwoCFluidState uses TypeTag

        // CompositionalFluidState
        Dumux::CompositionalFluidState<Scalar, FluidSystem> compositionalFluidState;
        collectedExceptions = Dumux::checkFluidState<Scalar>(compositionalFluidState);
        if (!collectedExceptions.empty())
        {
            std::cout<<"Dumux::CompositionalFluidState: \n"<<collectedExceptions<<"\n";
            success++;
        }

        // ImmiscibleFluidState
        Dumux::ImmiscibleFluidState<Scalar, FluidSystem> immiscibleFluidState;
        collectedExceptions = Dumux::checkFluidState<Scalar>(immiscibleFluidState);
        if (!collectedExceptions.empty())
        {
            std::cout<<"Dumux::ImmiscibleFluidState: \n"<<collectedExceptions<<"\n";
            success++;
        }

        // IsothermalImmiscibleFluidState
        Dumux::IsothermalImmiscibleFluidState<Scalar, FluidSystem> isothermalImmiscibleFluidState;
        collectedExceptions = Dumux::checkFluidState<Scalar>(isothermalImmiscibleFluidState);
        if (!collectedExceptions.empty())
        {
            std::cout<<"Dumux::IsothermalImmiscibleFluidState: \n"<<collectedExceptions<<"\n";
            success++;
        }

        // NonEquilibriumFluidState
        Dumux::NonEquilibriumFluidState<Scalar, FluidSystem> nonEquilibriumFluidState;
        collectedExceptions = Dumux::checkFluidState<Scalar>(nonEquilibriumFluidState);
        if (!collectedExceptions.empty())
        {
            std::cout<<"Dumux::NonEquilibriumFluidState: \n"<<collectedExceptions<<"\n";
            success++;
        }

        // NonEquilibriumEnergyFluidState uses TypeTag

        // NonEquilibriumMassFluidState uses TypeTag

        // PressureOverlayFluidState
        Dumux::PressureOverlayFluidState<Scalar, BaseFluidState> pressureOverlayFluidState(baseFs);
        collectedExceptions = Dumux::checkFluidState<Scalar>(pressureOverlayFluidState);
        if (!collectedExceptions.empty())
        {
            std::cout<<"Dumux::PressureOverlayFluidState: \n"<<collectedExceptions<<"\n";
            success++;
        }

        // SaturationOverlayFluidState
        Dumux::SaturationOverlayFluidState<Scalar, BaseFluidState> saturationOverlayFluidState(baseFs);
        collectedExceptions = Dumux::checkFluidState<Scalar>(saturationOverlayFluidState);
        if (!collectedExceptions.empty())
        {
            std::cout<<"Dumux::SaturationOverlayFluidState: \n"<<collectedExceptions<<"\n";
            success++;
        }

        // TemperatureOverlayFluidState
        Dumux::TemperatureOverlayFluidState<Scalar, BaseFluidState> temperatureOverlayFluidState(baseFs);
        collectedExceptions = Dumux::checkFluidState<Scalar>(temperatureOverlayFluidState);
        if (!collectedExceptions.empty())
        {
            std::cout<<"Dumux::TemperatureOverlayFluidState: \n"<<collectedExceptions<<"\n";
            success++;
        }
    }


    //////////////////////////
    // check all fluid systems

    // 1p
    {   typedef Dumux::FluidSystems::OneP<Scalar, Liquid> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   typedef Dumux::FluidSystems::OneP<Scalar, Gas> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // 2p-immiscible
    {   typedef Dumux::FluidSystems::TwoPImmiscible<Scalar, Liquid, Liquid> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   typedef Dumux::FluidSystems::TwoPImmiscible<Scalar, Liquid, Gas> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   typedef Dumux::FluidSystems::TwoPImmiscible<Scalar, Gas, Liquid> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // base

    // Brine -- Air
    {   const bool enableComplexRelations=false;
        typedef Dumux::FluidSystems::BrineAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   const bool enableComplexRelations=true;
        typedef Dumux::FluidSystems::BrineAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // Brine -- CO2

    // H2O -- Air
    {   typedef Dumux::SimpleH2O<Scalar> H2O;
        const bool enableComplexRelations=false;
        typedef Dumux::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   typedef Dumux::SimpleH2O<Scalar> H2O;
        const bool enableComplexRelations=true;
        typedef Dumux::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   typedef Dumux::H2O<Scalar> H2O;
        const bool enableComplexRelations=false;
        typedef Dumux::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   typedef Dumux::H2O<Scalar> H2O;
        const bool enableComplexRelations=true;
        typedef Dumux::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // gas phase
    {   typedef Dumux::GasPhase<Scalar, H2O> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   typedef Dumux::FluidSystems::GasPhase<Scalar, H2O> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Mesitylene
    {   typedef Dumux::FluidSystems::H2OAirMesitylene<Scalar> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Xylene
    {   typedef Dumux::FluidSystems::H2OAirXylene<Scalar> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2
    {   typedef Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   typedef Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/true> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2 -- kinetic
    {   typedef Dumux::FluidSystems::H2ON2Kinetic<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   typedef Dumux::FluidSystems::H2ON2Kinetic<Scalar, /*enableComplexRelations=*/true> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2 -- o2
    {   typedef Dumux::FluidSystems::H2ON2O2<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   typedef Dumux::FluidSystems::H2ON2O2<Scalar, /*enableComplexRelations=*/true> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // liquid phase
    {   typedef Dumux::FluidSystems::LiquidPhase<Scalar, H2O> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   typedef Dumux::FluidSystems::LiquidPhase<Scalar, H2O> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // pure water simple
    {   typedef Dumux::FluidSystems::PureWaterSimpleFluidSystem<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }
    {   typedef Dumux::FluidSystems::PureWaterSimpleFluidSystem<Scalar, /*enableComplexRelations=*/true> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // spe5
    {   typedef Dumux::FluidSystems::Spe5<Scalar> FluidSystem;
        success += Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    return success;
}
