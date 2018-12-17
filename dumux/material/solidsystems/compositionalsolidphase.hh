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
 * \ingroup SolidSystems
 * \brief A solid phase consisting of multiple inert solid components.
 */
#ifndef DUMUX_SOLIDSYSTEMS_COMPOSITIONAL_SOLID_PHASE_HH
#define DUMUX_SOLIDSYSTEMS_COMPOSITIONAL_SOLID_PHASE_HH

#include <string>
#include <dune/common/exceptions.hh>

namespace Dumux {
namespace SolidSystems {

/*!
 * \ingroup SolidSystems
 * \brief A solid phase consisting of multiple inert solid components.
 * \note a solid is considered inert if it cannot dissolve in a liquid and
 *       and cannot increase its mass by precipitation from a fluid phase.
 * \note inert components have to come after all non-inert components
 */
template <class Scalar, class Component1, class Component2, int numInert = 0>
class CompositionalSolidPhase
{
public:
    using ComponentOne = Component1;
    using ComponentTwo = Component2;


    /****************************************
     * Solid phase related static parameters
     ****************************************/
    static constexpr int numComponents = 2;
    static constexpr int numInertComponents = numInert;
    static constexpr int comp0Idx = 0;
    static constexpr int comp1Idx = 1;


    /*!
     * \brief Return the human readable name of a solid phase
     *
     * \param compIdx The index of the solid phase to consider
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx)
        {
            case comp0Idx: return ComponentOne::name();
            case comp1Idx: return ComponentTwo::name();
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief A human readable name for the solid system.
     */
    static std::string name()
    { return "s"; }

    /*!
     * \brief Returns whether the phase is incompressible
     */
    static constexpr bool isCompressible(int compIdx)
    { return false; }

    /*!
     * \brief Returns whether the component is inert (doesn't react)
     */
    static constexpr bool isInert()
    { return (numComponents == numInertComponents); }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the component.
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx)
        {
            case comp0Idx: return ComponentOne::molarMass();
            case comp1Idx: return ComponentTwo::molarMass();
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the solid phase at a given pressure and temperature.
     */
    template <class SolidState>
    static Scalar density(const SolidState& solidState)
    {
        Scalar rho1 = ComponentOne::solidDensity(solidState.temperature());
        Scalar rho2 = ComponentTwo::solidDensity(solidState.temperature());
        Scalar volFrac1 = solidState.volumeFraction(comp0Idx);
        Scalar volFrac2 = solidState.volumeFraction(comp1Idx);

        return (rho1*volFrac1+
               rho2*volFrac2)/(volFrac1+volFrac2);
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the solid phase at a given pressure and temperature.
     */
    template <class SolidState>
    static Scalar density(const SolidState& solidState, const int compIdx)
    {
        switch (compIdx)
        {
            case comp0Idx: return ComponentOne::solidDensity(solidState.temperature());
            case comp1Idx: return ComponentTwo::solidDensity(solidState.temperature());
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief The molar density of the solid phase at a given pressure and temperature.
     */
    template <class SolidState>
    static Scalar molarDensity(const SolidState& solidState, const int compIdx)
    {
        switch (compIdx)
        {
            case comp0Idx: return ComponentOne::solidDensity(solidState.temperature())/ComponentOne::molarMass();
            case comp1Idx: return ComponentTwo::solidDensity(solidState.temperature())/ComponentTwo::molarMass();
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief Thermal conductivity of the solid \f$\mathrm{[W/(m K)]}\f$.
     */
    template <class SolidState>
    static Scalar thermalConductivity(const SolidState &solidState)
    {
        Scalar lambda1 = ComponentOne::solidThermalConductivity(solidState.temperature());
        Scalar lambda2 = ComponentTwo::solidThermalConductivity(solidState.temperature());
        Scalar volFrac1 = solidState.volumeFraction(comp0Idx);
        Scalar volFrac2 = solidState.volumeFraction(comp1Idx);

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
        Scalar volFrac1 = solidState.volumeFraction(comp0Idx);
        Scalar volFrac2 = solidState.volumeFraction(comp1Idx);

        return (c1*volFrac1+
               c2*volFrac2)/(volFrac1+volFrac2);
    }

};

} // end namespace SolidSystems
} // end namespace Dumux

#endif
