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
 * \ingroup Fluidsystems
 * \brief @copybrief Dumux::FluidSystems::Base
 */
#ifndef DUMUX_BASE_FLUID_SYSTEM_HH
#define DUMUX_BASE_FLUID_SYSTEM_HH

#include <string>

#include <dune/common/exceptions.hh>
#include <dumux/common/typetraits/typetraits.hh>
#include "nullparametercache.hh"

namespace Dumux {
namespace FluidSystems {

/*!
* \ingroup Fluidsystems
* \brief Fluid system base class.
*
* \note Always derive your fluid system from this class to be sure
*       that all basic functionality is available!
*/
template <class ScalarType, class Implementation>
class Base
{
public:
    //! export the scalar type
    using Scalar = ScalarType;

    //! The type of parameter cache objects
    using ParameterCache = NullParameterCache;

    /*!
     * \brief Some properties of the fluid system
     */
    // \{

    //! If the fluid system only contains tracer components
    static constexpr bool isTracerFluidSystem()
    { return false; }

    /*!
     * \brief Get the main component of a given phase if possible
     *
     * \param phaseIdx The index of the fluid phase to consider
     * \todo Unfortunately we currently still have the assumption in some volume variables (e.g. 1pnc, 2pnc)
     *       that the main component index of a phase is equal to the phase index of that phase. This means
     *       changing this only works if the volume variables are written accordingly.
     * \note This only makes sense if this is not a tracer fluid system (then the bulk component is not balanced)
     */
    template<class I = Implementation, std::enable_if_t<!I::isTracerFluidSystem(), int> = 0>
    static constexpr int getMainComponent(int phaseIdx)
    { return phaseIdx; }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be compressible.
     *
     * Compressible means that the partial derivative of the density
     * to the fluid pressure is always larger than zero.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    template<class T = Implementation>
    static constexpr bool isCompressible(int phaseIdx)
    {
        static_assert(AlwaysFalse<T>::value, "Mandatory function not implemented: isCompressible(phaseIdx)");
        return true;
    }

    /*!
     * \brief Returns whether the fluids are miscible
     */
    template<class T = Implementation>
    static constexpr bool isMiscible()
    {
        static_assert(AlwaysFalse<T>::value, "Mandatory function not implemented: isMiscible()");
        return true;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        have a constant viscosity.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool viscosityIsConstant(int phaseIdx)
    { return false; }

    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    { return "DefaultPhaseName"; }

    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string componentName(int phaseIdx)
    { return "DefaultComponentName"; }

    // \}

    /*!
     * \brief Calculate the density \f$\mathrm{[kg/m^3]}\f$ of a fluid phase
     * \param fluidState The fluid state
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "density() method not implemented by the fluid system.");
    }

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
     * \brief Calculate the molar density \f$\mathrm{[mol/m^3]}\f$ of a fluid phase
     * \param fluidState The fluid state
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState,
                          int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "molarDensity() method not implemented by the fluid system.");
    }

    /*!
     * \brief Calculate the molar density \f$\mathrm{[mol/m^3]}\f$ of a fluid phase
     * \param fluidState The fluid state
     * \param paramCache mutable parameters
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        return Implementation::molarDensity(fluidState, phaseIdx);
    }

    /*!
     * \brief Calculate the fugacity coefficient \f$\mathrm{[Pa]}\f$ of an individual
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
     * \param phaseIdx Index of the fluid phase
     * \param compIdx Index of the component
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "fugacityCoefficient() method not implemented by the fluid system.");
    }

    /*!
     * \brief Calculate the fugacity coefficient \f$\mathrm{[Pa]}\f$ of an individual
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
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "viscosity() method not implemented by the fluid system.");
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
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "diffusionCoefficient() method not implemented by the fluid system.");
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
     * \param phaseIdx Index of the fluid phase
     * \param compIIdx Index of the component i
     * \param compJIdx Index of the component j
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)

    {
        DUNE_THROW(Dune::NotImplemented, "binaryDiffusionCoefficient() method not implemented by the fluid system.");
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
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                                 int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "enthalpy() method not implemented by the fluid system.");
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
     * \param phaseIdx Index of the fluid phase
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "thermalConductivity() method not implemented by the fluid system.");
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
                               int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "heatCapacity() method not implemented by the fluid system.");
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

} // end namespace FluidSystems

} // end namespace Dumux

#endif
