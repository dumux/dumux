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
 * \brief This file provides the actual code for the fluid systems test.
 *
 * It is not directly in test_fluidsystems.cc so that external modules
 * like dumux-devel can use it easily.
 */

#ifndef DUMUX_CHECK_FLUIDSYSTEM_HH
#define DUMUX_CHECK_FLUIDSYSTEM_HH

#include <exception>
#include <string>

#include <dune/common/classname.hh>

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
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/fluidstates/isothermalimmiscible.hh>
#include <dumux/material/fluidstates/nonequilibrium.hh>
#include <dumux/material/fluidstates/nonequilibriummass.hh>
#include <dumux/material/fluidstates/pressureoverlay.hh>
#include <dumux/material/fluidstates/pseudo1p2c.hh>
#include <dumux/material/fluidstates/saturationoverlay.hh>
#include <dumux/material/fluidstates/temperatureoverlay.hh>

namespace Dumux
{

/*!
 * \brief This fluid state ensures that only the allowed quantities are accessed.
 */
template<class ScalarType, class FluidSystem, class BaseFluidState = CompositionalFluidState<ScalarType, FluidSystem> >
class HairSplittingFluidState: protected BaseFluidState
{
public:
    //! Export the type used for scalars
    using typename BaseFluidState::Scalar;

    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;

    HairSplittingFluidState()
    {
        // set some fake values
        BaseFluidState::setTemperature(293.15);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            BaseFluidState::setSaturation(phaseIdx, 1.0 / numPhases);
            BaseFluidState::setPressure(phaseIdx, 1e5);
            BaseFluidState::setDensity(phaseIdx, 1.0);
            BaseFluidState::setMolarDensity(phaseIdx, 1.0);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                BaseFluidState::setMoleFraction(phaseIdx, compIdx, 1.0 / numComponents);
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

    Scalar wettingPhase() const
    {
        assert(allowComposition_);
        return BaseFluidState::wettingPhase();
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
            std::cout << "HairSplittingFluidState: pressure called but not allowed" << std::endl;
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
int checkFluidState(const BaseFluidState &fs)
{
    std::cout << "Testing fluid state '" << Dune::className<BaseFluidState>() << "'\n";

    // fluid states must be copy-able
    BaseFluidState tmpFs(fs);
    tmpFs = fs;

    // output strings
    std::string collectedErrors;
    std::string collectedWarnings;

    // make sure the fluid state provides all mandatory methods
    Scalar DUNE_UNUSED val;

    try
    {
        val = fs.temperature(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedErrors += "error: fluidState.temperature() throws exception!\n";
    }
    try
    {
        val = fs.pressure(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedErrors += "error: fluidState.pressure() throws exception!\n";
    }
    try
    {
        val = fs.moleFraction(/*phaseIdx=*/0, /*compIdx=*/0);
    } catch (...)
    {
        collectedErrors += "error: fluidState.moleFraction() throws exception!\n";
    }
    try
    {
        val = fs.massFraction(/*phaseIdx=*/0, /*compIdx=*/0);
    } catch (...)
    {
        collectedErrors += "error: fluidState.massFraction() throws exception!\n";
    }
    try
    {
        val = fs.averageMolarMass(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedErrors += "error: fluidState.averageMolarMass() throws exception!\n";
    }
    try
    {
        val = fs.molarity(/*phaseIdx=*/0, /*compIdx=*/0);
    } catch (...)
    {
        collectedErrors += "error: fluidState.molarity() throws exception!\n";
    }
    try
    {
        val = fs.molarDensity(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedErrors += "error: fluidState.molarDensity() throws exception!\n";
    }
    try
    {
        val = fs.molarVolume(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedErrors += "error: fluidState.molarVolume() throws exception!\n";
    }
    try
    {
        val = fs.density(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedErrors += "error: fluidState.density() throws exception!\n";
    }
    try
    {
        val = fs.saturation(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedErrors += "error: fluidState.saturation() throws exception!\n";
    }
    try
    {
        val = fs.fugacity(/*phaseIdx=*/0, /*compIdx=*/0);
    } catch (...)
    {
        collectedErrors += "error: fluidState.fugacity() throws exception!\n";
    }
    try
    {
        val = fs.fugacityCoefficient(/*phaseIdx=*/0, /*compIdx=*/0);
    } catch (...)
    {
        collectedErrors += "error: fluidState.fugacityCoefficient() throws exception!\n";
    }
    try
    {
        val = fs.enthalpy(/*phaseIdx=*/0);
    } catch (Dune::NotImplemented&)
    {
        collectedWarnings += "warning: fluidState.enthalpy() is not implemented\n";
    } catch (...)
    {
        collectedErrors += "error: fluidState.enthalpy() throws exception!\n";
    }
    try
    {
        val = fs.internalEnergy(/*phaseIdx=*/0);
    } catch (Dune::NotImplemented&)
    {
        collectedWarnings += "warning: fluidState.internalEnergy() is not implemented\n";
    } catch (...)
    {
        collectedErrors += "error: fluidState.internalEnergy() throws exception!\n";
    }
    try
    {
        val = fs.viscosity(/*phaseIdx=*/0);
    } catch (...)
    {
        collectedErrors += "error: fluidState.viscosity() throws exception!\n";
    }

    std::cout << collectedErrors;
//     std::cout << collectedWarnings;
    if (collectedErrors.empty()) // success
    {
        std::cout << "... successfull" << std::endl;
        std::cout << "----------------------------------" << std::endl;
        return 0;
    }
    else
    {
        std::cout << "... failed" << std::endl;
        std::cout << "----------------------------------" << std::endl;
        return 1;
    }
}

/*!
 * \brief This is a consistency check for FluidSystems.
 *
 * \param enablePhaseRestriction Parameter passed to the fluidState. If set to true,
 *        the fluidState will only allow calls to properties of the current phase.
 * \note While this is very common, it is not necessarily the case for all FluidSystems.
 *       We keep this, because it might help finding mistakes in FluidSystems that have this invariant.
 *       If you verified that a fluid system does not have this invariant you can set this option to false.
 */
template<class Scalar, class FluidSystem>
int checkFluidSystem(bool enablePhaseRestriction = true)
{
    int success = 0;
    std::cout << "Testing fluid system '" << Dune::className<FluidSystem>() << "'\n";
    FluidSystem::init();

    // output strings
    std::string collectedErrors;
    std::string collectedWarnings;

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
    using PC = typename FluidSystem::ParameterCache;
    PC paramCache;
    try
    {
        paramCache.updateAll(fs);
    } catch (...)
    {
        collectedErrors += "error: paramCache.updateAll() throws exception!\n";
    }
    try
    {
        paramCache.updateAll(fs, /*except=*/PC::None);
    } catch (...)
    {
        collectedErrors += "error: paramCache.updateAll(none) throws exception!\n";
    }
    try
    {
        paramCache.updateAll(fs, /*except=*/PC::Temperature | PC::Pressure | PC::Composition);
    } catch (...)
    {
        collectedErrors += "error: paramCache.updateAll(T, p, x) throws exception!\n";
    }
    try
    {
        paramCache.updateAllPressures(fs);
    } catch (...)
    {
        collectedErrors += "error: paramCache.updateAllPressures() throws exception!\n";
    }

    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
    {
        if (enablePhaseRestriction)
            fs.restrictToPhase(phaseIdx);
        try
        {
            paramCache.updatePhase(fs, phaseIdx);
        } catch (...)
        {
            collectedErrors += "error: paramCache.updatePhase() throws exception!\n";
        }
        try
        {
            paramCache.updatePhase(fs, phaseIdx, /*except=*/PC::None);
        } catch (...)
        {
            collectedErrors += "error: paramCache.updatePhase(none) throws exception!\n";
        }
        try
        {
            paramCache.updatePhase(fs, phaseIdx, /*except=*/PC::Temperature | PC::Pressure | PC::Composition);
        } catch (...)
        {
            collectedErrors += "error: paramCache.updatePhase(T, p , x) throws exception!\n";
        }
        try
        {
            paramCache.updateTemperature(fs, phaseIdx);
        } catch (...)
        {
            collectedErrors += "error: paramCache.updateTemperature() throws exception!\n";
        }
        try
        {
            paramCache.updatePressure(fs, phaseIdx);
        } catch (...)
        {
            collectedErrors += "error: paramCache.updatePressure() throws exception!\n";
        }
        try
        {
            paramCache.updateComposition(fs, phaseIdx);
        } catch (...)
        {
            collectedErrors += "error: paramCache.updateComposition() throws exception!\n";
        }
        try
        {
            paramCache.updateSingleMoleFraction(fs, phaseIdx, /*compIdx=*/0);
        } catch (...)
        {
            collectedErrors += "error: paramCache.updateSingleMoleFraction() throws exception!\n";
        }
    }

    // some value to make sure the return values of the fluid system
    // are convertible to scalars
    Scalar DUNE_UNUSED val;

    // actually check the fluid system API
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
    {
        fs.allowPressure(FluidSystem::isCompressible(phaseIdx));
        fs.allowComposition(true);
        fs.allowDensity(false);
        if (enablePhaseRestriction)
            fs.restrictToPhase(phaseIdx);
        try
        {
            val = FluidSystem::density(fs, paramCache, phaseIdx);
        } catch (const std::exception& e)
        {
            collectedErrors += "error: FluidSystem::density() throws exception: " + std::string(e.what()) + "\n";
        }
        try
        {
            val = FluidSystem::molarDensity(fs, paramCache, phaseIdx);
        } catch (const std::exception& e)
        {
            collectedErrors += "error: FluidSystem::molarDensity() throws exception: " + std::string(e.what()) + "\n";
        }
        fs.allowPressure(true);
        fs.allowDensity(true);
        try
        {
            val = FluidSystem::viscosity(fs, paramCache, phaseIdx);
        } catch (const std::exception& e)
        {
            collectedErrors += "error: FluidSystem::viscosity() throws exception: " + std::string(e.what()) + "\n";
        }
        try
        {
            val = FluidSystem::enthalpy(fs, paramCache, phaseIdx);
        } catch (Dune::NotImplemented&)
        {
            collectedWarnings += "warning: FluidSystem::enthalpy() is not implemented\n";
        } catch (const std::exception& e)
        {
            collectedErrors += "error: FluidSystem::enthalpy() throws exception: " + std::string(e.what()) + "\n";
        }
        try
        {
            val = FluidSystem::heatCapacity(fs, paramCache, phaseIdx);
        } catch (Dune::NotImplemented&)
        {
            collectedWarnings += "warning: FluidSystem::heatCapacity() is not implemented\n";
        } catch (const std::exception& e)
        {
            collectedErrors += "error: FluidSystem::heatCapacity() throws exception: " + std::string(e.what()) + "\n";
        }
        try
        {
            val = FluidSystem::thermalConductivity(fs, paramCache, phaseIdx);
        } catch (Dune::NotImplemented&)
        {
            collectedWarnings += "warning: FluidSystem::thermalConductivity() is not implemented\n";
        } catch (const std::exception& e)
        {
            collectedErrors += "error: FluidSystem::thermalConductivity() throws exception: " + std::string(e.what()) + "\n";
        }

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            fs.allowComposition(!FluidSystem::isIdealMixture(phaseIdx));
            try
            {
                val = FluidSystem::fugacityCoefficient(fs, paramCache, phaseIdx, compIdx);
            } catch (Dune::NotImplemented&)
            {
                collectedWarnings += "warning: FluidSystem::fugacityCoefficient() is not implemented\n";
            } catch (const std::exception& e)
            {
                collectedErrors += "error: FluidSystem::fugacityCoefficient() throws exception: " + std::string(e.what()) + "\n";
            }
            fs.allowComposition(true);
            try
            {
                val = FluidSystem::diffusionCoefficient(fs, paramCache, phaseIdx, compIdx);
            } catch (Dune::NotImplemented&)
            {
                collectedWarnings += "warning: FluidSystem::diffusionCoefficient() is not implemented\n";
            } catch (Dune::InvalidStateException&)
            {
                collectedWarnings += "warning: FluidSystem::diffusionCoefficient() gives invalid state exception\n";
            } catch (const std::exception& e)
            {
                collectedErrors += "error: FluidSystem::diffusionCoefficient() throws exception: " + std::string(e.what()) + "\n";
            }
            for (int comp2Idx = 0; comp2Idx < numComponents; ++comp2Idx)
            {
                try
                {
                    val = FluidSystem::binaryDiffusionCoefficient(fs, paramCache, phaseIdx, compIdx, comp2Idx);
                } catch (Dune::NotImplemented&)
                {
                    collectedWarnings += "warning: FluidSystem::binaryDiffusionCoefficient() is not implemented\n";
                } catch (Dune::InvalidStateException&)
                {
                    collectedWarnings += "warning: FluidSystem::binaryDiffusionCoefficient() gives invalid state exception\n";
                } catch (const std::exception& e)
                {
                    collectedErrors += "error: FluidSystem::binaryDiffusionCoefficient() throws exception: " + std::string(e.what()) + "\n";
                }
            }
        }
    }

    // test for phaseName(), isGas() and isIdealGas()
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
    {
        std::string
        DUNE_UNUSED name = FluidSystem::phaseName(phaseIdx);
        bool DUNE_UNUSED
        bVal = FluidSystem::isGas(phaseIdx);
        bVal = FluidSystem::isIdealGas(phaseIdx);
    }

    // test for componentName()
    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
    {
        val = FluidSystem::molarMass(compIdx);
        std::string
        DUNE_UNUSED name = FluidSystem::componentName(compIdx);
    }

    std::cout << collectedErrors;
//     std::cout << collectedWarnings;
    if (collectedErrors.empty()) // success
    {
        std::cout << "... successful" << std::endl;
        std::cout << "----------------------------------" << std::endl;
        return 0;
    }
    else
    {
        std::cout << "... failed" << std::endl;
        std::cout << "----------------------------------" << std::endl;
        return 1;
    }
    std::cout << "----------------------------------\n";
    return success;
}

} // end namespace Dumux

#endif
