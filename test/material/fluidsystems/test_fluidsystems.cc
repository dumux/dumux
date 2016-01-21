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
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/h2on2liquidphase.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidsystems/h2oairmesitylene.hh>
#include <dumux/material/fluidsystems/h2oairxylene.hh>

// include all fluid states
#include <dumux/material/fluidstates/pressureoverlay.hh>
#include <dumux/material/fluidstates/saturationoverlay.hh>
#include <dumux/material/fluidstates/temperatureoverlay.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/fluidstates/nonequilibrium.hh>
#include <dumux/material/fluidstates/immiscible.hh>

int main()
{
    typedef double Scalar;
    typedef Dumux::H2O<Scalar> H2O;
    typedef Dumux::N2<Scalar> N2;

    typedef Dumux::LiquidPhase<Scalar, H2O> Liquid;
    typedef Dumux::GasPhase<Scalar, N2> Gas;

    // check all fluid states
    {
        typedef Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false> FluidSystem;

        // CompositionalFluidState
        {   Dumux::CompositionalFluidState<Scalar, FluidSystem> fs;
            std::string collectedExceptions = Dumux::checkFluidState<Scalar>(fs);
            if (!collectedExceptions.empty()){
                std::cout<<"Dumux::CompositionalFluidState: \n"<<collectedExceptions<<"\n";
            }}

        // NonEquilibriumFluidState
        {   Dumux::NonEquilibriumFluidState<Scalar, FluidSystem> fs;
            std::string collectedExceptions = Dumux::checkFluidState<Scalar>(fs);
            if (!collectedExceptions.empty()){
                std::cout<<"Dumux::NonEquilibriumFluidState: \n"<<collectedExceptions<<"\n";
            }}

        // ImmiscibleFluidState
        {   Dumux::ImmiscibleFluidState<Scalar, FluidSystem> fs;
            std::string collectedExceptions = Dumux::checkFluidState<Scalar>(fs);
            if (!collectedExceptions.empty()){
                std::cout<<"Dumux::ImmiscibleFluidState: \n"<<collectedExceptions<<"\n";
            }}

        // IsothermalImmiscibleFluidState
        {   Dumux::IsothermalImmiscibleFluidState<Scalar, FluidSystem> fs;
            std::string collectedExceptions = Dumux::checkFluidState<Scalar>(fs);
            if (!collectedExceptions.empty()){
                std::cout<<"Dumux::IsothermalImmiscibleFluidState: \n"<<collectedExceptions<<"\n";
            }}

        typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> BaseFluidState;
        BaseFluidState baseFs;

        // TemperatureOverlayFluidState
        {   Dumux::TemperatureOverlayFluidState<Scalar, BaseFluidState> fs(baseFs);
            std::string collectedExceptions = Dumux::checkFluidState<Scalar>(fs);
            if (!collectedExceptions.empty()){
                std::cout<<"Dumux::TemperatureOverlayFluidState: \n"<<collectedExceptions<<"\n";
            } }

        // PressureOverlayFluidState
        {   Dumux::PressureOverlayFluidState<Scalar, BaseFluidState> fs(baseFs);
            std::string collectedExceptions = Dumux::checkFluidState<Scalar>(fs);
            if (!collectedExceptions.empty()){
                std::cout<<"Dumux::PressureOverlayFluidState: \n"<<collectedExceptions<<"\n";
            }}

        // SaturationOverlayFluidState
        {   Dumux::SaturationOverlayFluidState<Scalar, BaseFluidState> fs(baseFs);
            std::string collectedExceptions = Dumux::checkFluidState<Scalar>(fs);
            if (!collectedExceptions.empty()){
                std::cout<<"Dumux::SaturationOverlayFluidState: \n"<<collectedExceptions<<"\n";
            }}
    }

    // H2O -- N2
    {   typedef Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/true> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2 -- liquid phase
    {   typedef Dumux::FluidSystems::H2ON2LiquidPhase<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::FluidSystems::H2ON2LiquidPhase<Scalar, /*enableComplexRelations=*/true> FluidSystem;
         Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air
    {   typedef Dumux::SimpleH2O<Scalar> H2O;
        const bool enableComplexRelations=false;
        typedef Dumux::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::SimpleH2O<Scalar> H2O;
        const bool enableComplexRelations=true;
        typedef Dumux::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::H2O<Scalar> H2O;
        const bool enableComplexRelations=false;
        typedef Dumux::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::H2O<Scalar> H2O;
        const bool enableComplexRelations=true;
        typedef Dumux::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Mesitylene
    {   typedef Dumux::FluidSystems::H2OAirMesitylene<Scalar> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Xylene
    {   typedef Dumux::FluidSystems::H2OAirXylene<Scalar> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // 2p-immiscible
    {   typedef Dumux::FluidSystems::TwoPImmiscible<Scalar, Liquid, Liquid> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::FluidSystems::TwoPImmiscible<Scalar, Liquid, Gas> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    {  typedef Dumux::FluidSystems::TwoPImmiscible<Scalar, Gas, Liquid> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    // 1p
    {   typedef Dumux::FluidSystems::OneP<Scalar, Liquid> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::FluidSystems::OneP<Scalar, Gas> FluidSystem;
        Dumux::checkFluidSystem<Scalar, FluidSystem>(); }

    return 0;
}
