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
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system without using
 *        any assumptions.
 */
#ifndef DUMUX_NONEQUILIBRIUM_FLUID_STATE_HH
#define DUMUX_NONEQUILIBRIUM_FLUID_STATE_HH

#include <cmath>
#include <algorithm>

#include <dune/common/exceptions.hh>

#include <dumux/common/valgrind.hh>

namespace Dumux
{

/*!
 * \ingroup FluidStates
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system without using
 *        any assumptions.
 */
template <class Scalar, class FluidSystem>
class NonEquilibriumFluidState
{
public:
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;

    NonEquilibriumFluidState()
    {
        // set the composition to 0
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)  {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                moleFraction_[phaseIdx][compIdx] = 0;

            averageMolarMass_[phaseIdx] = 0;
            sumMoleFractions_[phaseIdx] = 0;
        }

        // make everything undefined so that valgrind will complain
        Valgrind::SetUndefined(*this);
    }

    /*****************************************************
     * Generic access to fluid properties (No assumptions
     * on thermodynamic equilibrium required)
     *****************************************************/
    /*!
     * @copydoc Dumux::CompositionalFluidState::saturation()
     */
    Scalar saturation(int phaseIdx) const
    { return saturation_[phaseIdx]; }

    /*!
     * @copydoc Dumux::CompositionalFluidState::moleFraction()
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    { return moleFraction_[phaseIdx][compIdx]; }


    /*!
     * @copydoc Dumux::CompositionalFluidState::massFraction()
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    {
      using std::abs;
      using std::max;
        return
            abs(sumMoleFractions_[phaseIdx])
            * moleFraction_[phaseIdx][compIdx]
            * FluidSystem::molarMass(compIdx)
            / max(1e-40, abs(averageMolarMass_[phaseIdx]));
    }

    /*!
     * @copydoc Dumux::CompositionalFluidState::averageMolarMass()
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return averageMolarMass_[phaseIdx]; }

    /*!
     * @copydoc Dumux::CompositionalFluidState::molarity()
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    { return molarDensity(phaseIdx)*moleFraction(phaseIdx, compIdx); }

    /*!
     * @copydoc Dumux::CompositionalFluidState::fugacityCoefficient()
     */
    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    { return fugacityCoefficient_[phaseIdx][compIdx]; }

    /*!
     * @copydoc Dumux::CompositionalFluidState::fugacity()
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    {
            return pressure_[phaseIdx]*fugacityCoefficient_[phaseIdx][compIdx]*moleFraction_[phaseIdx][compIdx];
    }

    Scalar fugacity(int compIdx) const
    {
           return fugacity(0, compIdx);
    }

    /*!
     * @copydoc Dumux::CompositionalFluidState::molarVolume()
     */
    Scalar molarVolume(int phaseIdx) const
    { return 1/molarDensity(phaseIdx); }

    /*!
     * @copydoc Dumux::CompositionalFluidState::density()
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * @copydoc Dumux::CompositionalFluidState::molarDensity()
     */
    Scalar molarDensity(int phaseIdx) const
    { return density_[phaseIdx]/averageMolarMass(phaseIdx); }


    /*!
     * @copydoc Dumux::CompositionalFluidState::temperature()
     */
    Scalar temperature(const int phaseIdx) const
    {  return temperature_[phaseIdx];  }

    /*!
     * \brief Get the equilibrium temperature \f$\mathrm{[K]}\f$ of the fluid phases.
     */
    Scalar temperature() const
    {
        return temperatureEquil_ ;
    }

    /*!
     * @copydoc Dumux::CompositionalFluidState::pressure()
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; }
    /*!
     * @copydoc Dumux::CompositionalFluidState::enthalpy()
     */
    Scalar enthalpy(int phaseIdx) const
    { return enthalpy_[phaseIdx]; }

    /*!
     * @copydoc Dumux::CompositionalFluidState::internalEnergy()
     */
    Scalar internalEnergy(int phaseIdx) const
    {
        return enthalpy_[phaseIdx]
                         - pressure(phaseIdx)/density(phaseIdx);
    }

    /*!
     * @copydoc Dumux::CompositionalFluidState::viscosity()
     */
    Scalar viscosity(int phaseIdx) const
    { return viscosity_[phaseIdx]; }

    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/
    /*!
     * \brief Set the temperature \f$\mathrm{[K]}\f$ of a fluid phase
     */
    void setTemperature(int phaseIdx, Scalar value)
    {
        temperature_[phaseIdx] = value;
    }

    /*!
     * \brief Set the temperature \f$\mathrm{[K]}\f$ of all fluid phases.
     */
    void setTemperature(Scalar value)
    {
        temperatureEquil_ = value;
    }

    /*!
     * \brief Set the fluid pressure of a phase \f$\mathrm{[Pa]}\f$
     */
    void setPressure(int phaseIdx, Scalar value)
    { pressure_[phaseIdx] = value; }

    /*!
     * \brief Set the saturation of a phase \f$\mathrm{[-]}\f$
     */
    void setSaturation(int phaseIdx, Scalar value)
    { saturation_[phaseIdx] = value; }

    /*!
     * \brief Set the mole fraction of a component in a phase \f$\mathrm{[-]}\f$
     */
    void setMoleFraction(int phaseIdx, int compIdx, Scalar value)
    {
        Valgrind::CheckDefined(value);
        Valgrind::SetDefined(sumMoleFractions_[phaseIdx]);
        Valgrind::SetDefined(averageMolarMass_[phaseIdx]);
        Valgrind::SetDefined(moleFraction_[phaseIdx][compIdx]);

        using std::isfinite;
        if (isfinite(averageMolarMass_[phaseIdx])) {
            Scalar delta = value - moleFraction_[phaseIdx][compIdx];

            moleFraction_[phaseIdx][compIdx] = value;

            sumMoleFractions_[phaseIdx] += delta;
            averageMolarMass_[phaseIdx] += delta*FluidSystem::molarMass(compIdx);
        }
        else {
            moleFraction_[phaseIdx][compIdx] = value;

            // re-calculate the mean molar mass
            sumMoleFractions_[phaseIdx] = 0.0;
            averageMolarMass_[phaseIdx] = 0.0;
            for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
                sumMoleFractions_[phaseIdx] += moleFraction_[phaseIdx][compJIdx];
                averageMolarMass_[phaseIdx] += moleFraction_[phaseIdx][compJIdx]*FluidSystem::molarMass(compJIdx);
            }
        }
    }

    /*!
     * \brief Set the fugacity of a component in a phase \f$\mathrm{[-]}\f$
     */
    void setFugacityCoefficient(int phaseIdx, int compIdx, Scalar value)
    { fugacityCoefficient_[phaseIdx][compIdx] = value; }

    /*!
     * \brief Set the density of a phase \f$\mathrm{[kg / m^3]}\f$
     */
    void setDensity(int phaseIdx, Scalar value)
    { density_[phaseIdx] = value; }

    /*!
     * \brief Set the specific enthalpy of a phase \f$\mathrm{[J/m^3]}\f$
     */
    void setEnthalpy(int phaseIdx, Scalar value)
    { enthalpy_[phaseIdx] = value; }

    /*!
     * \brief Set the dynamic viscosity of a phase \f$\mathrm{[Pa s]}\f$
     */
    void setViscosity(int phaseIdx, Scalar value)
    { viscosity_[phaseIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     * \param fs Fluidstate
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            averageMolarMass_[phaseIdx] = 0;
            sumMoleFractions_[phaseIdx] = 0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                moleFraction_[phaseIdx][compIdx] = fs.moleFraction(phaseIdx, compIdx);
                fugacityCoefficient_[phaseIdx][compIdx] = fs.fugacityCoefficient(phaseIdx, compIdx);
                averageMolarMass_[phaseIdx] += moleFraction_[phaseIdx][compIdx]*FluidSystem::molarMass(compIdx);
                sumMoleFractions_[phaseIdx] += moleFraction_[phaseIdx][compIdx];
            }

            pressure_[phaseIdx] = fs.pressure(phaseIdx);
            saturation_[phaseIdx] = fs.saturation(phaseIdx);
            density_[phaseIdx] = fs.density(phaseIdx);
            enthalpy_[phaseIdx] = fs.enthalpy(phaseIdx);
            viscosity_[phaseIdx] = fs.viscosity(phaseIdx);
            temperature_[phaseIdx] = fs.temperature(phaseIdx);
        }
    }

    /*!
     * \brief Make sure that all attributes are defined.
     *
     * This method does not do anything if the program is not run
     * under valgrind. If it is, then valgrind will print an error
     * message if some attributes of the object have not been properly
     * defined.
     */
    void checkDefined() const
    {
#if HAVE_VALGRIND && ! defined NDEBUG
        for (int i = 0; i < numPhases; ++i) {
            for (int j = 0; j < numComponents; ++j) {
                Valgrind::CheckDefined(fugacityCoefficient_[i][j]);
                Valgrind::CheckDefined(moleFraction_[i][j]);
            }
            Valgrind::CheckDefined(averageMolarMass_[i]);
            Valgrind::CheckDefined(pressure_[i]);
            Valgrind::CheckDefined(saturation_[i]);
            Valgrind::CheckDefined(density_[i]);
            Valgrind::CheckDefined(temperature_[i]);
            Valgrind::CheckDefined(enthalpy_[i]);
            Valgrind::CheckDefined(viscosity_[i]);
        }
#endif // HAVE_VALGRIND
    }

protected:
    Scalar moleFraction_[numPhases][numComponents];
    Scalar fugacityCoefficient_[numPhases][numComponents];

    Scalar averageMolarMass_[numPhases];
    Scalar sumMoleFractions_[numPhases];
    Scalar pressure_[numPhases];
    Scalar saturation_[numPhases];
    Scalar density_[numPhases];
    Scalar enthalpy_[numPhases];
    Scalar viscosity_[numPhases];
    Scalar temperature_[numPhases]; //
    Scalar temperatureEquil_;
};

} // end namespace Dumux

#endif
