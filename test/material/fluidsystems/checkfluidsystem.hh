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
 * \brief This file provides the actual code for the fluid systems
 *        test.
 *
 * It is not directly in test_fluidsystems.cc so that external modules
 * like dumux-devel can use it easily
 */
#ifndef DUMUX_CHECK_FLUIDSYSTEM_HH
#define DUMUX_CHECK_FLUIDSYSTEM_HH

#include <dune/common/classname.hh>

// include all fluid systems in dumux-stable
#include <dumux/material/fluidsystems/1pfluidsystem.hh>
#include <dumux/material/fluidsystems/2pimmisciblefluidsystem.hh>
#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <dumux/material/fluidsystems/h2on2liquidphasefluidsystem.hh>
#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>
#include <dumux/material/fluidsystems/h2oairmesitylenefluidsystem.hh>
#include <dumux/material/fluidsystems/h2oairxylenefluidsystem.hh>
#include <dumux/material/fluidsystems/spe5fluidsystem.hh>

// include all fluid states
#include <dumux/material/fluidstates/pressureoverlayfluidstate.hh>
#include <dumux/material/fluidstates/saturationoverlayfluidstate.hh>
#include <dumux/material/fluidstates/temperatureoverlayfluidstate.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>
#include <dumux/material/fluidstates/nonequilibriumfluidstate.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include <dumux/material/fluidstates/isothermalimmisciblefluidstate.hh>

namespace Dumux
{

/*! \brief This fluid state ensures that only the allowed quantities
 * are accessed
 */
template<class Scalar, class FluidSystem, class BaseFluidState = Dumux::CompositionalFluidState<Scalar, FluidSystem> >
class HairSplittingFluidState: protected BaseFluidState
{
public:
    enum
    {
        numPhases = FluidSystem::numPhases
    };
    enum
    {
        numComponents = FluidSystem::numComponents
    };

    HairSplittingFluidState()
    {
        // set some fake values
        BaseFluidState::setTemperature(293.15);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            BaseFluidState::setSaturation(phaseIdx, 1.0 / numPhases);
            BaseFluidState::setDensity(phaseIdx, 1.0);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                BaseFluidState::setMoleFraction(phaseIdx, compIdx, 1.0 / numComponents);

            }
        }

        // initially, do not allow anything
        allowTemperature(false);
        allowPressure(false);
        allowComposition(false);
        allowDensity(false);

        // do not allow accessing any phase
        restrictToPhase(1000);
    }

    void allowTemperature(bool yesno)
    {
        allowTemperature_ = yesno;
    }

    void allowPressure(bool yesno)
    {
        allowPressure_ = yesno;
    }

    void allowComposition(bool yesno)
    {
        allowComposition_ = yesno;
    }

    void allowDensity(bool yesno)
    {
        allowDensity_ = yesno;
    }

    void restrictToPhase(int phaseIdx)
    {
        restrictPhaseIdx_ = phaseIdx;
    }

    Scalar temperature(int phaseIdx) const
    {
        assert(allowTemperature_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        return BaseFluidState::temperature(phaseIdx);
    }

    Scalar partialPressure(int phaseIdx, int compIdx) const
    {
        assert(allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        return BaseFluidState::partialPressure(phaseIdx, compIdx);
    }

    Scalar pressure(int phaseIdx) const
    {
        if (!allowPressure_)
        {
            DUNE_THROW(Dune::InvalidStateException,
                       "pressure called but not allowed\n\n");
        }
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        return BaseFluidState::pressure(phaseIdx);
    }

    Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        assert(allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        return BaseFluidState::moleFraction(phaseIdx, compIdx);
    }

    Scalar massFraction(int phaseIdx, int compIdx) const
    {
        assert(allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        return BaseFluidState::massFraction(phaseIdx, compIdx);
    }

    Scalar averageMolarMass(int phaseIdx) const
    {
        assert(allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        return BaseFluidState::averageMolarMass(phaseIdx);
    }

    Scalar molarity(int phaseIdx, int compIdx) const
    {
        assert(allowDensity_ && allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        return BaseFluidState::molarity(phaseIdx, compIdx);
    }

    Scalar molarDensity(int phaseIdx) const
    {
        assert(allowDensity_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        return BaseFluidState::molarDensity(phaseIdx);
    }

    Scalar molarVolume(int phaseIdx) const
    {
        assert(allowDensity_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        return BaseFluidState::molarVolume(phaseIdx);
    }

    Scalar density(int phaseIdx) const
    {
        assert(allowDensity_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        return BaseFluidState::density(phaseIdx);
    }

    Scalar saturation(int phaseIdx) const
    {
        assert(false);
        return BaseFluidState::saturation(phaseIdx);
    }

    Scalar fugacity(int phaseIdx, int compIdx) const
    {
        assert(false);
        return BaseFluidState::fugacity(phaseIdx, compIdx);
    }

    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    {
        assert(false);
        return BaseFluidState::fugacityCoefficient(phaseIdx, compIdx);
    }

    Scalar enthalpy(int phaseIdx) const
    {
        assert(false);
        return BaseFluidState::enthalpy(phaseIdx);
    }

    Scalar internalEnergy(int phaseIdx) const
    {
        assert(false);
        return BaseFluidState::internalEnergy(phaseIdx);
    }

    Scalar viscosity(int phaseIdx) const
    {
        assert(false);
        return BaseFluidState::viscosity(phaseIdx);
    }

private:
    bool allowSaturation_;
    bool allowTemperature_;
    bool allowPressure_;
    bool allowComposition_;
    bool allowDensity_;
    int restrictPhaseIdx_;
};

template<class Scalar, class BaseFluidState>
std::string checkFluidState(const BaseFluidState &fs)
{
    // fluid states must be copy-able
    BaseFluidState tmpFs(fs);
    tmpFs = fs;

    // a fluid state must provide a checkDefined() method
    fs.checkDefined();

    std::string collectedExceptions;
    // make sure the fluid state provides all mandatory methods
    Scalar DUNE_UNUSED val;

    try
    {
        val = fs.temperature(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.temperature() throws exception!\n";
    }
    try
    {
        val = fs.pressure(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.pressure() throws exception!\n";
    }
    try
    {
        val = fs.moleFraction(/*phaseIdx=*/0, /*compIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.moleFraction() throws exception!\n";
    }
    try
    {
        val = fs.massFraction(/*phaseIdx=*/0, /*compIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.massFraction() throws exception!\n";
    }
    try
    {
        val = fs.averageMolarMass(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.averageMolarMass() throws exception!\n";
    }
    try
    {
        val = fs.molarity(/*phaseIdx=*/0, /*compIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.molarity() throws exception!\n";
    }
    try
    {
        val = fs.molarDensity(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.molarDensity() throws exception!\n";
    }
    try
    {
        val = fs.molarVolume(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.molarVolume() throws exception!\n";
    }
    try
    {
        val = fs.density(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.density() throws exception!\n";
    }
    try
    {
        val = fs.saturation(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.saturation() throws exception!\n";
    }
    try
    {
        val = fs.fugacity(/*phaseIdx=*/0, /*compIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.fugacity() throws exception!\n";
    }
    try
    {
        val = fs.fugacityCoefficient(/*phaseIdx=*/0, /*compIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.fugacityCoefficient() throws exception!\n";
    }
    try
    {
        val = fs.enthalpy(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.enthalpy() throws exception!\n";
    }
    try
    {
        val = fs.internalEnergy(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.internalEnergy() throws exception!\n";
    }
    try
    {
        val = fs.viscosity(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedExceptions += "fluidState.viscosity() throws exception!\n";
    }
    return collectedExceptions;
}

template<class Scalar, class FluidSystem>
void checkFluidSystem()
{
    std::cout << "Testing fluid system '" << Dune::className<FluidSystem>() << "'\n";

    // make sure the fluid system provides the number of phases and
    // the number of components
    enum
    {
        numPhases = FluidSystem::numPhases
    };
    enum
    {
        numComponents = FluidSystem::numComponents
    };

    HairSplittingFluidState<Scalar, FluidSystem> fs;
    fs.allowTemperature(true);
    fs.allowPressure(true);
    fs.allowComposition(true);
    fs.restrictToPhase(-1);

    // check whether the parameter cache adheres to the API
    typedef typename FluidSystem::ParameterCache PC;
    PC paramCache;
    try
    {
        paramCache.updateAll(fs);
    } catch (...)
    {
    }
    try
    {
        paramCache.updateAll(fs, /*except=*/PC::None);
    } catch (...)
    {
    }
    try
    {
        paramCache.updateAll(fs, /*except=*/PC::Temperature | PC::Pressure | PC::Composition);
    } catch (...)
    {
    }
    try
    {
        paramCache.updateAllPressures(fs);
    } catch (...)
    {
    };

    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
    {
        fs.restrictToPhase(phaseIdx);
        try
        {
            paramCache.updatePhase(fs, phaseIdx);
        } catch (...)
        {
        }
        try
        {
            paramCache.updatePhase(fs, phaseIdx, /*except=*/PC::None);
        } catch (...)
        {
        }
        try
        {
            paramCache.updatePhase(fs, phaseIdx, /*except=*/PC::Temperature | PC::Pressure | PC::Composition);
        } catch (...)
        {
        }
        try
        {
            paramCache.updateTemperature(fs, phaseIdx);
        } catch (...)
        {
        }
        try
        {
            paramCache.updatePressure(fs, phaseIdx);
        } catch (...)
        {
        }
        try
        {
            paramCache.updateComposition(fs, phaseIdx);
        } catch (...)
        {
        }
        try
        {
            paramCache.updateSingleMoleFraction(fs, phaseIdx, /*compIdx=*/0);
        } catch (...)
        {
        }
    }

    // some value to make sure the return values of the fluid system
    // are convertible to scalars
    Scalar DUNE_UNUSED val;

    // actually check the fluid system API
    FluidSystem::init();
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
    {
        fs.restrictToPhase(phaseIdx);
        fs.allowPressure(FluidSystem::isCompressible(phaseIdx));
        fs.allowComposition(true);
        fs.allowDensity(false);
        try
        {
            val = FluidSystem::density(fs, paramCache, phaseIdx);
        } catch (Dune::Exception e)
        {
            std::cout << "\ndensity calculation throws exception:\n" << e.what();
        }

        fs.allowPressure(true);
        fs.allowDensity(true);
        try
        {
            val = FluidSystem::viscosity(fs, paramCache, phaseIdx);
        } catch (...)
        {
        }
        try
        {
            val = FluidSystem::enthalpy(fs, paramCache, phaseIdx);
        } catch (...)
        {
        }
        try
        {
            val = FluidSystem::heatCapacity(fs, paramCache, phaseIdx);
        } catch (...)
        {
        }
        try
        {
            val = FluidSystem::thermalConductivity(fs, paramCache, phaseIdx);
        } catch (...)
        {
        }

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            fs.allowComposition(!FluidSystem::isIdealMixture(phaseIdx));
            try
            {
                val = FluidSystem::fugacityCoefficient(fs, paramCache, phaseIdx, compIdx);
            } catch (...)
            {
            }
            fs.allowComposition(true);
            try
            {
                val = FluidSystem::diffusionCoefficient(fs, paramCache, phaseIdx, compIdx);
            } catch (...)
            {
            }
            for (int comp2Idx = 0; comp2Idx < numComponents; ++comp2Idx)
            {
                try
                {
                    val = FluidSystem::binaryDiffusionCoefficient(fs, paramCache, phaseIdx, compIdx, comp2Idx);
                } catch (...)
                {
                }
            }
        }
    }

    // test for phaseName(), isLiquid() and isIdealGas()
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
    {
        std::string
        DUNE_UNUSED name = FluidSystem::phaseName(phaseIdx);
        bool DUNE_UNUSED
        bVal = FluidSystem::isLiquid(phaseIdx);
        bVal = FluidSystem::isIdealGas(phaseIdx);
    }

    // test for componentName()
    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
    {
        val = FluidSystem::molarMass(compIdx);
        std::string
        DUNE_UNUSED name = FluidSystem::componentName(compIdx);
    }

    std::cout << "----------------------------------\n";
}

} // end namespace Dumux

#endif // DUMUX_CHECK_FLUIDSYSTEM_HH
