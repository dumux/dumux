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
 * \brief @copybrief Dumux::FluidSystems::BaseFluidSystem
 *
 */
#ifndef DUMUX_BASE_FLUID_SYSTEM_HH
#define DUMUX_BASE_FLUID_SYSTEM_HH

#include "nullparametercache.hh"

namespace Dumux
{
namespace FluidSystems
{
/*!
 * \ingroup Fluidsystems
 * \brief Fluid system base class.
 *
 * \note This fluid system neglects the contribution of gas-molecules on the
 *       liquid phase. This contribution is probably not big, but somebody
 *       has to find out its influence.
 */
template <class Scalar, class Implementation>
class BaseFluidSystem
{
public:
    //! The type of parameter cache objects
    typedef NullParameterCache ParameterCache;

    /*!
     * \brief Calculate the density \f$\mathrm{[kg/m^3]}\f$ of a fluid phase
     * \param fluidState The fluid state
     * \param paramCache mutable parameters
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        return Implementation::density(fluidState, phaseIdx);
    }

    /*!
     * \brief Calculate the fugacity coefficient \f$\mathrm{[-]}\f$ of an individual
     *        component in a fluid phase
     *
     * The fugacity coefficient \f$\mathrm{\phi^\kappa_\alpha}\f$ is connected to the
     * fugacity \f$\mathrm{f^\kappa_\alpha}\f$ and the component's mole
     * fraction \f$\mathrm{x^\kappa_\alpha}\f$ by means of the relation
     *
     * \f[
     f^\kappa_\alpha = \phi^\kappa_\alpha\;x^\kappa_\alpha\;p_\alpha
     * \f]
     *
     * \param fluidState The fluid state
     * \param paramCache mutable parameters
     * \param phaseIdx Index of the fluid phase
     * \param compIdx Index of the component
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx,
                                      int compIdx)
    {
        return Implementation::fugacityCoefficient(fluidState, phaseIdx, compIdx);
    }

    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
     * \param fluidState The fluid state
     * \param paramCache mutable parameters
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const ParameterCache &paramCache,
                            int phaseIdx)
    {
        return Implementation::viscosity(fluidState, phaseIdx);
    }

    /*!
     * \brief Calculate the binary molecular diffusion coefficient for
     *        a component in a fluid phase \f$\mathrm{[mol^2 * s / (kg*m^3)]}\f$
     * \param fluidState The fluid state
     * \param paramCache mutable parameters
     * \param phaseIdx Index of the fluid phase
     * \param compIdx Index of the component
     * Molecular diffusion of a component \f$\mathrm{\kappa}\f$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D \mathbf{grad} \mu_\kappa \f]
     *
     * where \f$\mathrm{\mu_\kappa}\f$ is the component's chemical potential,
     * \f$\mathrm{D}\f$ is the diffusion coefficient and \f$\mathrm{J}\f$ is the
     * diffusive flux. \f$\mathrm{\mu_\kappa}\f$ is connected to the component's
     * fugacity \f$\mathrm{f_\kappa}\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$\mathrm{p_\alpha}\f$ and \f$\mathrm{T_\alpha}\f$ are the fluid phase'
     * pressure and temperature.
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       const ParameterCache &paramCache,
                                       int phaseIdx,
                                       int compIdx)
    {
        return Implementation::diffusionCoefficient(fluidState, phaseIdx, compIdx);
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
     *        \f$\mathrm{i}\f$ and \f$\mathrm{j}\f$ in this phase.
     * \param fluidState The fluid state
     * \param paramCache mutable parameters
     * \param phaseIdx Index of the fluid phase
     * \param compIIdx Index of the component i
     * \param compJIdx Index of the component j
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             const ParameterCache &paramCache,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)

    {
        return Implementation::binaryDiffusionCoefficient(fluidState, phaseIdx, compIIdx, compJIdx);
    }

    /*!
     * \brief Given a phase's composition, temperature, pressure and
     *        density, calculate its specific enthalpy \f$\mathrm{[J/kg]}\f$.
     * \param fluidState The fluid state
     * \param paramCache mutable parameters
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                                 const ParameterCache &paramCache,
                                 int phaseIdx)
    {
        return Implementation::enthalpy(fluidState, phaseIdx);
    }

    /*!
     * \brief Thermal conductivity \f$\lambda_\alpha \f$ of a fluid phase \f$\mathrm{[W/(m K)]}\f$.
     * \param fluidState The fluid state
     * \param paramCache mutable parameters
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx)
    {
        return Implementation::thermalConductivity(fluidState, phaseIdx);
    }

    /*!
     * \brief Specific isobaric heat capacity \f$c_{p,\alpha}\f$ of a fluid phase \f$\mathrm{[J/(kg*K)]}\f$.
     *
     * \param fluidState represents all relevant thermodynamic quantities of a fluid system
     * \param paramCache mutable parameters
     * \param phaseIdx Index of the fluid phase
     *
     * Given a fluid state, an up-to-date parameter cache and a phase index, this method
     * computes the isobaric heat capacity \f$c_{p,\alpha}\f$ of the fluid phase. The isobaric
     * heat capacity is defined as the partial derivative of the specific enthalpy \f$h_\alpha\f$
     * to the fluid pressure \f$p_\alpha\f$:
     *
     * \f$ c_{p,\alpha} = \frac{\partial h_\alpha}{\partial p_\alpha} \f$
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               const ParameterCache &paramCache,
                               int phaseIdx)
    {
        return Implementation::heatCapacity(fluidState, phaseIdx);
    }
};

} // end namespace
} // end namespace

#endif
