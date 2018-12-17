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
 * \ingroup SolidStates
 * \brief Represents all relevant thermodynamic quantities of a inert solid system.
 */
#ifndef DUMUX_INERT_SOLID_STATE_HH
#define DUMUX_INERT_SOLID_STATE_HH

namespace Dumux {

/*!
 * \ingroup SolidStates
 * \brief Represents all relevant thermodynamic quantities of a inert solid system.
 */
template <class Scalar, class SolidSystemType>
class InertSolidState
{
public:
    using SolidSystem = SolidSystemType;

    enum
    {
        numComponents = SolidSystem::numComponents,
        numInertComponents = SolidSystem::numInertComponents,
    };

    /*!
     * \brief Allows compile-time evaluation of if the solid system
     *        is inert or takes part in any kind of reactions.
     */
    static constexpr bool isInert()
    {
        static_assert(SolidSystem::isInert(), "Only inert solid systems are allowed with the InertSolidState");
        return true;
    }

    /*!
     * \brief The average molar mass \f$\overline M_\alpha\f$ of phase \f$\alpha\f$ in \f$\mathrm{[kg/mol]}\f$
     *
     * Since this is an inert InertSolidState we simply consider the molarMass of the
     * pure component/phase.
     */
    Scalar averageMolarMass() const
    { return SolidSystem::molarMass(); }

    /*!
     * \brief The porosity of the porous medium
     */
    Scalar porosity() const
    {
        Scalar sumVolumeFraction = 0.0;
        for (int compIdx =0; compIdx < numComponents; ++compIdx)
            sumVolumeFraction += volumeFraction(compIdx);
        Scalar porosity = 1-sumVolumeFraction;
        return porosity;
    }

    //! The mass density of the solid phase in \f$\mathrm{[kg/m^3]}\f$
    Scalar density() const { return density_; }

    //! The heat capacity of the solid phase in \f$\mathrm{[J/(kg*K)}\f$
    Scalar heatCapacity() const { return heatCapacity_; }

    //! The thermal conductivity of the solid phase in \f$\mathrm{[[W/(m*K)]}\f$
    Scalar thermalConductivity() const { return thermalConducivity_; }

    //! The temperature of the solid phase in \f$\mathrm{[K]}\f$
    Scalar temperature() const { return temperature_; }

    //! The volume fraction of a solid component within the solid phase
    Scalar volumeFraction(const int compIdx) const { return volumeFraction_[compIdx]; }


    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a solid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the mass density \f$\rho_\alpha\f$ and the mean molar mass \f$\overline M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{\overline M_\alpha} \;.\f]
     */
    Scalar molarDensity() const
    { return density_/averageMolarMass(); }

   /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic InertSolidState interface but specific for each
     * implementation...
     *****************************************************/

    /*!
     * \brief Retrieve all parameters from an arbitrary solid state.
     * \param sst The inert solid state
     *
     * \note If the other solid state object is inconsistent with the
     *       thermodynamic equilibrium, the result of this method is
     *       undefined.
     */
    template <class SolidState>
    void assign(const SolidState &sst)
    {
        temperature_ = sst.temperature();
        density_ = sst.density();
        thermalConducivity_ = sst.thermalConductivity();
        heatCapacity_ = sst.heatCapacity();
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            volumeFraction_[compIdx] = sst.volumeFraction(compIdx);
    }

    /*!
     * \brief Set the temperature \f$\mathrm{[K]}\f$  of the solid phase
     */
    void setTemperature(Scalar value)
    { temperature_ = value; }

    /*!
     * \brief Set the density of the solid phase
     */
    void setDensity(Scalar value)
    { density_ = value; }

    /*!
     * \brief Set the heat capacity of the solid phase
     */
    void setThermalConductivity(Scalar value)
    { thermalConducivity_ = value; }

    /*!
     * \brief Set the thermal conductivity of the solid phase
     */
    void setHeatCapacity(Scalar value)
    { heatCapacity_ = value; }

    /*!
     * \brief Set the volume fraction of a solid component
     */
    void setVolumeFraction(const int compIdx, Scalar value)
    { volumeFraction_[compIdx] = value; }

protected:
    Scalar density_;
    Scalar temperature_;
    Scalar volumeFraction_[numComponents];
    Scalar heatCapacity_;
    Scalar thermalConducivity_;
};

} // end namespace Dumux

#endif
