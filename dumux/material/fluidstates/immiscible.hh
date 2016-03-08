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
 *        multi-phase fluid system assuming immiscibility and
 *        thermodynamic equilibrium.
 */
#ifndef DUMUX_IMMISCIBLE_FLUID_STATE_HH
#define DUMUX_IMMISCIBLE_FLUID_STATE_HH

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
class ImmiscibleFluidState
{
public:
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int numPhases = FluidSystem::numPhases;
    static_assert(numPhases == numComponents,
                  "The number of phases must be equal to the number of "
                  "components if immiscibility is assumed!");

    ImmiscibleFluidState()
    { Valgrind::SetUndefined(*this); }

    template <class FluidState>
    ImmiscibleFluidState(FluidState &fs)
    { assign(fs); }

    /*****************************************************
     * Generic access to fluid properties (No assumptions
     * on thermodynamic equilibrium required)
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
     * @copydoc CompositionalFluidState::fugacity()
     *
     * To avoid numerical issues with code that assumes miscibility,
     * we return a fugacity of 0 for components which do not mix with
     * the specified phase. (Actually it is undefined, but for finite
     * fugacity coefficients, the only way to get components
     * completely out of a phase is 0 to feed it zero fugacity.)
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    {
        if (phaseIdx == compIdx)
            return pressure(phaseIdx);
        else
            return 0;
    }

    /*!
     *  @copydoc CompositionalFluidState::fugacityCoefficient()
     *
     * Since we assume immiscibility, the fugacity coefficients for
     * the components which are not miscible with the phase is
     * infinite. Beware that this will very likely break your code if
     * you don't keep that in mind.
     */
    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    {
        if (phaseIdx == compIdx)
            return 1.0;
        else
            return std::numeric_limits<Scalar>::infinity();
    }

    /*!
     * @copydoc CompositionalFluidState::partialPressure()
     *
     * To avoid numerical issues with code that assumes miscibility,
     * we return a partial pressure of 0 for components which do not mix with
     * the specified phase. Actually it is undefined.
     */
    Scalar partialPressure(int phaseIdx, int compIdx) const
    {
        if (phaseIdx == compIdx)
            return pressure(phaseIdx);
        else
            return 0;
    }

    /*!
     * @copydoc CompositionalFluidState::molarVolume()
     *
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
     * @copydoc CompositionalFluidState::temperature()
     */
    Scalar temperature(int phaseIdx) const
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
    { return enthalpy_[phaseIdx]; }

    /*!
     * @copydoc CompositionalFluidState::internalEnergy()
     */
    Scalar internalEnergy(int phaseIdx) const
    { return enthalpy_[phaseIdx] - pressure(phaseIdx)/density(phaseIdx); }

    /*!
     * @copydoc CompositionalFluidState::viscosity()
     */
    Scalar viscosity(int phaseIdx) const
    { return viscosity_[phaseIdx]; }

    /*****************************************************
     * Access to fluid properties which only make sense
     * if assuming thermodynamic equilibrium
     *****************************************************/

    /*!
     * \brief The temperature within the domain \f$\mathrm{[K]}\f$
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     * \brief The fugacity of a component  \f$\mathrm{[Pa]}\f$
     *
     * This assumes chemical equilibrium.
     */
    Scalar fugacity(int compIdx) const
    { return fugacity(0, compIdx); }


    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     * \param fs Fluidstate
     *
     * \note If the other fluid state object is inconsistent with the
     *       thermodynamic equilibrium, the result of this method is
     *       undefined.
     */
    template <class FluidState>
    void assign(const FluidState &fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            pressure_[phaseIdx] = fs.pressure(phaseIdx);
            saturation_[phaseIdx] = fs.saturation(phaseIdx);
            density_[phaseIdx] = fs.density(phaseIdx);
            enthalpy_[phaseIdx] = fs.enthalpy(phaseIdx);
            viscosity_[phaseIdx] = fs.viscosity(phaseIdx);
        }
        temperature_ = fs.temperature(0);
    }

    /*!
     * \brief Set the temperature \f$\mathrm{[K]}\f$ of a fluid phase
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
     * \brief Set the specific enthalpy of a phase \f$\mathrm{[J/kg]}\f$
     */
    void setEnthalpy(int phaseIdx, Scalar value)
    { enthalpy_[phaseIdx] = value; }

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
            //for (int j = 0; j < numComponents; ++j) {
            //    Valgrind::CheckDefined(fugacityCoefficient_[i][j]);
            //}
            Valgrind::CheckDefined(pressure_[i]);
            Valgrind::CheckDefined(saturation_[i]);
            Valgrind::CheckDefined(density_[i]);
            //Valgrind::CheckDefined(internalEnergy_[i]);
            Valgrind::CheckDefined(viscosity_[i]);
        }

        Valgrind::CheckDefined(temperature_);
#endif // HAVE_VALGRIND
    }

protected:
    Scalar pressure_[numPhases];
    Scalar saturation_[numPhases];
    Scalar density_[numPhases];
    Scalar enthalpy_[numPhases];
    Scalar viscosity_[numPhases];
    Scalar temperature_;
};

} // end namespace Dumux

#endif
