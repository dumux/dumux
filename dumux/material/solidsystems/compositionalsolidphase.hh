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
#ifndef DUMUX_SOLIDSYSTEMS_COMPOSITIONAL_SOLID_PHASE_HH
#define DUMUX_SOLIDSYSTEMS_COMPOSITIONAL_SOLID_PHASE_HH

namespace Dumux {
namespace SolidSystems {

/*!
 * \ingroup SolidSystems
 * \brief A solid phase consisting of a single inert solid component
 * \note a solid is considered inert if it can't dissolve in a liquid and
 *       and can't increase its mass by precipitation from a fluid phase.
 */
template <class Scalar, class Component1, bool isInert1, class Component2, bool isInert2>
class CompositionalSolidPhase
{
public:
    using ComponentOne = Component1;
    using ComponentTwo = Component2;


    /****************************************
     * Solid phase related static parameters
     ****************************************/
    static constexpr int numComponents = 2;
    static constexpr int numInertComponents = isInert1 ? (isInert2 ? 2 : 1) : (isInert2 ? 1 : 0);
    static constexpr int componentOneIdx = 0;
    static constexpr int componentTwoIdx = 1;


    /*!
     * \brief Return the human readable name of a solid phase
     *
     * \param phaseIdx The index of the solid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        switch (phaseIdx) {
        case componentOneIdx: return ComponentOne::name();
        case componentTwoIdx: return ComponentTwo::name();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief A human readable name for the component.
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        DUNE_THROW(Dune::InvalidStateException, "ComponentName does not exist " << compIdx);
    }

    /*!
     * \brief A human readable name for the solid system.
     */
    static std::string name()
    { return "CompositionalSolidPhase"; }

    /*!
     * \brief Returns whether the phase is incompressible
     */
    static constexpr bool isCompressible(int phaseIdx)
    { return false; }

    /*!
     * \brief Returns whether the component is inert (doesn't react)
     */
    static constexpr bool isInert()
    {
        if (numComponents == numInertComponents)
            return true;
        else
            return false;
    }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the component.
     */
    static Scalar molarMass(int phaseIdx)
    {
        switch (phaseIdx) {
        case componentOneIdx: return ComponentOne::molarMass();
        case componentTwoIdx: return ComponentTwo::molarMass();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the solid phase at a given pressure and temperature.
     */
    template <class SolidState>
    static Scalar density(const SolidState& solidState)
    {
        Scalar rho1 = ComponentOne::solidDensity(solidState.temperature());
        Scalar rho2 = ComponentTwo::solidDensity(solidState.temperature());
        Scalar volFrac1 = solidState.volumeFraction(componentOneIdx);
        Scalar volFrac2 = solidState.volumeFraction(componentTwoIdx);

        return (rho1*volFrac1+
               rho2*volFrac2)/(volFrac1+volFrac2);
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the solid phase at a given pressure and temperature.
     */
    template <class SolidState>
    static Scalar density(const SolidState& solidState, const int phaseIdx)
    {
        switch (phaseIdx)
        {
            case componentOneIdx: return ComponentOne::solidDensity(solidState.temperature());
            case componentTwoIdx: return ComponentTwo::solidDensity(solidState.temperature());
        }
         DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief The molar density of the solid phase at a given pressure and temperature.
     */
    template <class SolidState>
    static Scalar molarDensity(const SolidState& solidState, const int phaseIdx)
    {
        switch (phaseIdx)
        {
            case componentOneIdx: return ComponentOne::solidDensity(solidState.temperature())/ComponentOne::molarMass();
            case componentTwoIdx: return ComponentTwo::solidDensity(solidState.temperature())/ComponentTwo::molarMass();
        }
         DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Thermal conductivity of the solid \f$\mathrm{[W/(m K)]}\f$.
     */
    template <class SolidState>
    static Scalar thermalConductivity(const SolidState &solidState)
    {
        Scalar lambda1 = ComponentOne::solidThermalConductivity(solidState.temperature());
        Scalar lambda2 = ComponentTwo::solidThermalConductivity(solidState.temperature());
        Scalar volFrac1 = solidState.volumeFraction(componentOneIdx);
        Scalar volFrac2 = solidState.volumeFraction(componentTwoIdx);

        return (lambda1*volFrac1+
               lambda2*volFrac2)/(volFrac1+volFrac2);
    }

    /*!
     * \brief Specific isobaric heat capacity of the pure solids \f$\mathrm{[J/(kg K)]}\f$.
     */
    template <class SolidState>
    static Scalar heatCapacity(const SolidState &solidState)
    {
        Scalar c1 = ComponentOne::solidHeatCapacity(solidState.temperature());
        Scalar c2 = ComponentTwo::solidHeatCapacity(solidState.temperature());
        Scalar volFrac1 = solidState.volumeFraction(componentOneIdx);
        Scalar volFrac2 = solidState.volumeFraction(componentTwoIdx);

        return (c1*volFrac1+
               c2*volFrac2)/(volFrac1+volFrac2);
    }

};

} // end namespace SolidSystems
} // end namespace Dumux

#endif
