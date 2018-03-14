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
 * \ingroup FluidStates
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
     * \brief Returns the saturation \f$S_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$.
     *
     * The saturation is defined as the pore space occupied by the fluid divided by the total pore space:
     *  \f[S_\alpha := \frac{\phi \mathcal{V}_\alpha}{\phi \mathcal{V}}\f]
     *
     * \param phaseIdx the index of the phase
     */
    Scalar saturation(int phaseIdx) const
    { return saturation_[phaseIdx]; }

    /*!
     * \brief Returns the molar fraction \f$x^\kappa_\alpha\f$ of the component \f$\kappa\f$ in fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$.
     *
     * The molar fraction \f$x^\kappa_\alpha\f$ is defined as the ratio of the number of molecules
     * of component \f$\kappa\f$ and the total number of molecules of the phase \f$\alpha\f$.
     * They are set either 1 or 0 in a phase since this is an immiscible fluidstate.
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    { return (phaseIdx == compIdx)?1.0:0.0; }

    /*!
     * \brief Returns the mass fraction \f$X^\kappa_\alpha\f$ of component \f$\kappa\f$ in fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$.
     *
     * They are set either 1 or 0 in a phase since this is an immiscible fluidstate.
     *
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    { return (phaseIdx == compIdx)?1.0:0.0; }

    /*!
     * \brief The average molar mass \f$\overline M_\alpha\f$ of phase \f$\alpha\f$ in \f$\mathrm{[kg/mol]}\f$
     *
     * The average molar mass is the mean mass of a mole of the
     * fluid at current composition. It is defined as the sum of the
     * component's molar masses weighted by the current mole fraction:
     * \f[\mathrm{ \overline M_\alpha = \sum_\kappa M^\kappa x_\alpha^\kappa}\f]
     *
     * Since this is an immiscible fluidstate we simply consider the molarMass of the
     * pure component/phase.
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return FluidSystem::molarMass(/*compIdx=*/phaseIdx); }

    /*!
     * \brief The molar concentration \f$c^\kappa_\alpha\f$ of component \f$\kappa\f$ in fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * This quantity is usually called "molar concentration" or just
     * "concentration", but there are many other (though less common)
     * measures for concentration.
     *
     * http://en.wikipedia.org/wiki/Concentration
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    { return molarDensity(phaseIdx)*moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief The fugacity \f$f^\kappa_\alpha\f$ of component \f$\kappa\f$
     *  in fluid phase \f$\alpha\f$ in \f$\mathrm{[Pa]}\f$
     * @copydoc ImmiscibleFluidState::fugacity()
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
     * \brief The fugacity coefficient \f$\Phi^\kappa_\alpha\f$ of component \f$\kappa\f$ in fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$
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
     * \brief The molar volume \f$v_{mol,\alpha}\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[m^3/mol]}\f$
     *
     * This quantity is the inverse of the molar density.
     */
    Scalar molarVolume(int phaseIdx) const
    { return 1/molarDensity(phaseIdx); }

    /*!
     * \brief The mass density \f$\rho_\alpha\f$ of the fluid phase
     *  \f$\alpha\f$ in \f$\mathrm{[kg/m^3]}\f$
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the mass density \f$\rho_\alpha\f$ and the mean molar mass \f$\overline M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{\overline M_\alpha} \;.\f]
     */
    Scalar molarDensity(int phaseIdx) const
    { return molarDensity_[phaseIdx]; }

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
     * \brief The pressure \f$p_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[Pa]}\f$
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; }

    /*!
     * \brief The specific enthalpy \f$h_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[J/kg]}\f$
     * This is not defined for an isothermal fluidstate.
     */
    Scalar enthalpy(int phaseIdx) const
    { DUNE_THROW(Dune::NotImplemented,"No enthalpy() function defined for isothermal systems!"); }

    /*!
     * \brief The specific internal energy \f$u_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[J/kg]}\f$
     *
     * The specific internal energy is defined by the relation:
     * \f[u_\alpha = h_\alpha - \frac{p_\alpha}{\rho_\alpha}\f]
     * This is not defined for an isothermal fluidstate.
     */
    Scalar internalEnergy(int phaseIdx) const
    { DUNE_THROW(Dune::NotImplemented,"No internalEnergy() function defined for isothermal systems!"); }

    /*!
     * \brief The dynamic viscosity \f$\mu_\alpha\f$ of fluid phase \f$\alpha\f$ in \f$\mathrm{[Pa s]}\f$
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
            molarDensity_[phaseIdx] = fs.molarDensity(phaseIdx);
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
     * \brief Set the molar density of a phase \f$\mathrm{[kg/m^3]}\f$
     */
    void setMolarDensity(int phaseIdx, Scalar value)
    { molarDensity_[phaseIdx] = value; }

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
            Valgrind::CheckDefined(molarDensity_[i]);
            Valgrind::CheckDefined(viscosity_[i]);
        }

        Valgrind::CheckDefined(temperature_);
#endif // HAVE_VALGRIND
    }

protected:
    Scalar pressure_[numPhases];
    Scalar saturation_[numPhases];
    Scalar density_[numPhases];
    Scalar molarDensity_[numPhases];
    Scalar viscosity_[numPhases];
    Scalar temperature_;
};

} // end namespace Dumux

#endif
