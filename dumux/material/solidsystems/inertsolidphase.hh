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
 * \ingroup SolidSystems
 * \brief @copybrief Dumux::SolidSystems::InertSolidPhase
 */
#ifndef DUMUX_SOLIDSYSTEMS_INERT_SOLID_PHASE_HH
#define DUMUX_SOLIDSYSTEMS_INERT_SOLID_PHASE_HH

namespace Dumux {
namespace SolidSystems {

/*!
 * \ingroup SolidSystems
 * \brief A solid phase consisting of a single inert solid component
 * \note a solid is considered inert if it can't dissolve in a liquid and
 *       and can't increase its mass by precipitation from a fluid phase.
 */
template <class Scalar, class ComponentT>
class InertSolidPhase
{
public:
    using Component = ComponentT;

    /****************************************
     * Solid phase related static parameters
     ****************************************/
    static constexpr int numComponents = 1;
    static constexpr int numInertComponents = 1;

    enum { lastInertCompIdx = numInertComponents-1 };

    /*!
     * \brief Return the human readable name of a solid phase
     *
     * \param phaseIdx The index of the solid phase to consider
     */
    static std::string phaseName(int phaseIdx = 0)
    { return Component::name(); }

    /*!
     * \brief A human readable name for the component.
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx = 0)
    { return Component::name(); }

    /*!
     * \brief A human readable name for the solid system.
     */
    static std::string name()
    { return Component::name(); }

    /*!
     * \brief Returns whether the phase is incompressible
     */
    static constexpr bool isCompressible(int phaseIdx = 0)
    { return false; }

    /*!
     * \brief Returns whether the component is inert (doesn't react)
     */
    static constexpr bool isInert()
    {  return true; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the component.
     */
    static Scalar molarMass(int compIdx = 0)
    { return Component::molarMass(); }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the solid phase at a given pressure and temperature.
     */
    static Scalar density(Scalar temperature)
    { return Component::solidDensity(temperature); }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the solid phase at a given pressure and temperature.
     */
    template <class SolidState>
    static Scalar density(const SolidState& solidState)
    {
        return density(solidState.temperature());
    }

   /*!
     * \brief Thermal conductivity of the solid \f$\mathrm{[W/(m K)]}\f$.
     */
    static Scalar thermalConductivity(Scalar temperature)
    { return Component::solidThermalConductivity(temperature); }

    /*!
     * \brief Thermal conductivity of the solid \f$\mathrm{[W/(m K)]}\f$.
     */
    template <class SolidState>
    static Scalar thermalConductivity(const SolidState &solidState)
    {
        return thermalConductivity(solidState.temperature());
    }

    /*!
     * \brief Specific isobaric heat capacity of the solid \f$\mathrm{[J/(kg K)]}\f$.
     */
    static Scalar heatCapacity(Scalar temperature)
    { return Component::solidHeatCapacity(temperature); }

    /*!
     * \brief Specific isobaric heat capacity of the solid \f$\mathrm{[J/(kg K)]}\f$.
     */
    template <class SolidState>
    static Scalar heatCapacity(const SolidState &solidState)
    {
        return heatCapacity(solidState.temperature());
    }

};

} // end namespace SolidSystems
} // end namespace Dumux

#endif
