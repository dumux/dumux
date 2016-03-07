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
 * \brief Represents all relevant thermodynamic quantities of a isothermal immiscible
 *        multi-phase fluid system
 */
#ifndef DUMUX_ISOIMMISCIBLE_FLUID_STATE_HH
#define DUMUX_ISOIMMISCIBLE_FLUID_STATE_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/valgrind.hh>

#include <limits>

namespace Dumux
{
/*!
 * \ingroup FluidStates
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase fluid system assuming immiscibility and
 *        thermodynamic equilibrium.
 */
template <class Scalar, class FluidSystem>
class IsothermalImmiscibleFluidState
{
public:
    static constexpr int numPhases = FluidSystem::numPhases;

    IsothermalImmiscibleFluidState()
    { Valgrind::SetUndefined(*this); }

    template <class FluidState>
    IsothermalImmiscibleFluidState(FluidState &fs)
    { assign(fs); }

    /*****************************************************
     * Generic access to fluid properties
     *****************************************************/
    /*!
     * @copydoc CompositionalFluidState::saturation()
     */
    Scalar saturation(int phaseIdx) const
    { return saturation_[phaseIdx]; }

    /*!
     * @copydoc CompositionalFluidState::moleFraction()
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    { return (phaseIdx == compIdx)?1.0:0.0; }

    /*!
     * @copydoc CompositionalFluidState::massFraction()
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    { return (phaseIdx == compIdx)?1.0:0.0; }

    /*!
     * @copydoc CompositionalFluidState::averageMolarMass()
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return FluidSystem::molarMass(/*compIdx=*/phaseIdx); }

    /*!
     * @copydoc CompositionalFluidState::molarity()
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    { return molarDensity(phaseIdx)*moleFraction(phaseIdx, compIdx); }

    /*!
     * @copydoc ImmiscibleFluidState::fugacity()
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    {
        if (phaseIdx == compIdx)
            return pressure(phaseIdx);
        else
            return 0;
    }

    /*!
     * @copydoc ImmiscibleFluidState::fugacityCoefficient()
     */
    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    {
        if (phaseIdx == compIdx)
            return 1.0;
        else
            return std::numeric_limits<Scalar>::infinity();
    }

    /*!
     * @copydoc CompositionalFluidState::molarVolume()
     */
    Scalar molarVolume(int phaseIdx) const
    { return 1/molarDensity(phaseIdx); }

    /*!
     * @copydoc CompositionalFluidState::density()
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * @copydoc CompositionalFluidState::molarDensity()
     */
    Scalar molarDensity(int phaseIdx) const
    { return density_[phaseIdx]/averageMolarMass(phaseIdx); }

    /*!
     * \brief The temperature of a fluid phase \f$\mathrm{[K]}\f$
     */
    Scalar temperature(int phaseIdx) const
    { return temperature_; }

    /*!
     * \brief The temperature within the domain \f$\mathrm{[K]}\f$
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     *  @copydoc CompositionalFluidState::pressure()
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; }

    /*!
     * @copydoc CompositionalFluidState::enthalpy()
     */
    Scalar enthalpy(int phaseIdx) const
    { DUNE_THROW(Dune::NotImplemented,"No enthalpy() function defined for isothermal systems!"); }

    /*!
     * @copydoc CompositionalFluidState::internalEnergy()
     */
    Scalar internalEnergy(int phaseIdx) const
    { DUNE_THROW(Dune::NotImplemented,"No internalEnergy() function defined for isothermal systems!"); }

    /*!
     * @copydoc CompositionalFluidState::viscosity()
     */
    Scalar viscosity(int phaseIdx) const
    { return viscosity_[phaseIdx]; }

    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     * \param fs Fluidstate
     */
    template <class FluidState>
    void assign(const FluidState &fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            pressure_[phaseIdx] = fs.pressure(phaseIdx);
            saturation_[phaseIdx] = fs.saturation(phaseIdx);
            density_[phaseIdx] = fs.density(phaseIdx);
            viscosity_[phaseIdx] = fs.viscosity(phaseIdx);
        }
        temperature_ = fs.temperature(0);
    }

    /*!
     * \brief Set the temperature \f$\mathrm{[K]]}\f$ of a fluid phase
     */
    void setTemperature(Scalar value)
    { temperature_ = value; }

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
     * \brief Set the density of a phase \f$\mathrm{[kg/m^3]}\f$
     */
    void setDensity(int phaseIdx, Scalar value)
    { density_[phaseIdx] = value; }

    /*!
     * \brief Set the dynamic viscosity of a phase \f$\mathrm{[Pa s]}\f$
     */
    void setViscosity(int phaseIdx, Scalar value)
    { viscosity_[phaseIdx] = value; }

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
            Valgrind::CheckDefined(pressure_[i]);
            Valgrind::CheckDefined(saturation_[i]);
            Valgrind::CheckDefined(density_[i]);
            Valgrind::CheckDefined(viscosity_[i]);
        }

        Valgrind::CheckDefined(temperature_);
#endif // HAVE_VALGRIND
    }

protected:
    Scalar pressure_[numPhases];
    Scalar saturation_[numPhases];
    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];
    Scalar temperature_;
};

} // end namespace Dumux

#endif
