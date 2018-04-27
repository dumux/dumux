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
 * \ingroup SolidStates
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase fluid system assuming immiscibility and
 *        thermodynamic equilibrium.
 */
#ifndef DUMUX_SOLID_STATE_COMPOSITIONAL_HH
#define DUMUX_SOLID_STATE_COMPOSITIONAL_HH

#include <dumux/common/valgrind.hh>

#include <limits>
#include "basesolidstate.hh"

namespace Dumux
{

/*!
 * \ingroup SolidStates
 * \brief Represents all relevant thermodynamic quantities of a inert
 *        solid system
 */
template <class Scalar, class SolidSystem>
class CompositionalSolidState
{
public:
    enum
    {
        numComponents = SolidSystem::numComponents,
        numInertComponents = SolidSystem::numInertComponents,
    };
    static constexpr bool isInert() { return SolidSystem::isInert(); }
    /*!
     * \brief The average molar mass \f$\overline M_\alpha\f$ of phase \f$\alpha\f$ in \f$\mathrm{[kg/mol]}\f$
     *
     * Since this is an inert CompositionalSolidState we simply consider the molarMass of the
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
        return std::max(porosity, minPorosity_);
    }
    /*!
     * \brief The mass density \f$\rho_\alpha\f$ of the fluid phase
     *  \f$\alpha\f$ in \f$\mathrm{[kg/m^3]}\f$
     */
    Scalar density() const
    { return density_; }

    Scalar heatCapacity() const
    { return heatCapacity_; }

    Scalar thermalConductivity() const
    { return thermalConducivity_; }

    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the mass density \f$\rho_\alpha\f$ and the mean molar mass \f$\overline M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{\overline M_\alpha} \;.\f]
     */
    Scalar molarDensity() const
    { return density_/averageMolarMass(); }

    /*!
     * \brief The temperature within the domain \f$\mathrm{[K]}\f$
     */
    Scalar temperature() const
    { return temperature_; }

    Scalar volumeFraction(const int phaseIdx) const
    { return volumeFraction_[phaseIdx]; }

   /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic CompositionalSolidState interface but specific for each
     * implementation...
     *****************************************************/

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     * \param sst CompositionalSolidState
     *
     * \note If the other fluid state object is inconsistent with the
     *       thermodynamic equilibrium, the result of this method is
     *       undefined.
     */
    template <class CompositionalSolidState>
    void assign(const CompositionalSolidState &sst)
    {
        temperature_ = sst.temperature();
        density_ = sst.density();
    }

    /*!
     * \brief Set the temperature \f$\mathrm{[K]}\f$  of the solid phase
     */
    void setTemperature(Scalar value)
    { temperature_ = value; }

    /*!
     * \brief Set the density  of the solid phase
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

    void setVolumeFraction(const int compIdx, Scalar value)
    { volumeFraction_[compIdx] = value; }

    void setMinPorosity(Scalar value)
    { minPorosity_ = value; }

protected:
    Scalar density_;
    Scalar temperature_;
    Scalar volumeFraction_[numComponents];
    Scalar minPorosity_;
    Scalar heatCapacity_;
    Scalar thermalConducivity_;
};

} // end namespace Dumux

#endif
