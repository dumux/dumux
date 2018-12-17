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
 * \ingroup FluidStates
 * \brief This is a fluid state which allows to set the fluid
 *        temperatures and takes all other quantities from an other
 *        fluid state.
 */
#ifndef DUMUX_TEMPERATURE_OVERLAY_FLUID_STATE_HH
#define DUMUX_TEMPERATURE_OVERLAY_FLUID_STATE_HH

namespace Dumux {

/*!
 * \ingroup FluidStates
 * \brief This is a fluid state which allows to set the fluid
 *        temperatures and takes all other quantities from an other
 *        fluid state.
 */
template <class FluidState>
class TemperatureOverlayFluidState
{
public:
    static constexpr int numPhases = FluidState::numPhases;
    static constexpr int numComponents = FluidState::numComponents;

    //! export the scalar type
    using Scalar = typename FluidState::Scalar;

    /*!
     * \brief Constructor
     *
     * \param fs Fluidstate
     * The overlay fluid state copies the saturation from the argument,
     * so it initially behaves exactly like the underlying fluid
     * state.
     */
    TemperatureOverlayFluidState(const FluidState &fs)
    : fs_(&fs)
    {
        temperature_ = fs.temperature(/*phaseIdx=*/0);
    }

    TemperatureOverlayFluidState(Scalar T, const FluidState &fs)
    : fs_(&fs), temperature_(T)
    {}

    // copy & move constructor / assignment operators
    TemperatureOverlayFluidState(const TemperatureOverlayFluidState &fs) = default;
    TemperatureOverlayFluidState(TemperatureOverlayFluidState &&fs) = default;
    TemperatureOverlayFluidState& operator=(const TemperatureOverlayFluidState &fs) = default;
    TemperatureOverlayFluidState& operator=(TemperatureOverlayFluidState &&fs) = default;

    /*****************************************************
     * Generic access to fluid properties (No assumptions
     * on thermodynamic equilibrium required)
     *****************************************************/
    /*!
     *  @copydoc CompositionalFluidState::saturation()
     */
    Scalar saturation(int phaseIdx) const
    { return fs_->saturation(phaseIdx); }

    /*!
     *  @copydoc CompositionalFluidState::moleFraction()
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    { return fs_->moleFraction(phaseIdx, compIdx); }

    /*!
     * @copydoc CompositionalFluidState::massFraction()
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    { return fs_->massFraction(phaseIdx, compIdx); }

    /*!
     *  @copydoc CompositionalFluidState::averageMolarMass()
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return fs_->averageMolarMass(phaseIdx); }

    /*!
     *  @copydoc CompositionalFluidState::molarity()
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    { return fs_->molarity(phaseIdx, compIdx); }

    /*!
     *  @copydoc CompositionalFluidState::fugacity()
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    { return fs_->fugacity(phaseIdx, compIdx); }

    /*!
     *  @copydoc CompositionalFluidState::fugacityCoefficient()
     */
    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    { return fs_->fugacityCoefficient(phaseIdx, compIdx); }

    /*!
     *  @copydoc CompositionalFluidState::molarVolume()
     */
    Scalar molarVolume(int phaseIdx) const
    { return fs_->molarVolume(phaseIdx); }

    /*!
     *  @copydoc CompositionalFluidState::density()
     */
    Scalar density(int phaseIdx) const
    { return fs_->density(phaseIdx); }

    /*!
     *  @copydoc CompositionalFluidState::molarDensity()
     */
    Scalar molarDensity(int phaseIdx) const
    { return fs_->molarDensity(phaseIdx); }

    /*!
     * \brief The temperature of a fluid phase \f$\mathrm{[K]}\f$
     */
    Scalar temperature(int phaseIdx) const
    { return temperature_; }

    /*!
     *  @copydoc CompositionalFluidState::pressure()
     */
    Scalar pressure(int phaseIdx) const
    { return fs_->pressure(phaseIdx); }

    /*!
     *  @copydoc CompositionalFluidState::enthalpy()
     */
    Scalar enthalpy(int phaseIdx) const
    { return fs_->enthalpy(phaseIdx); }

    /*!
     *  @copydoc CompositionalFluidState::internalEnergy()
     */
    Scalar internalEnergy(int phaseIdx) const
    { return fs_->internalEnergy(phaseIdx); }

    /*!
     *  @copydoc CompositionalFluidState::viscosity()
     */
    Scalar viscosity(int phaseIdx) const
    { return fs_->viscosity(phaseIdx); }


    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/
    /*!
     * \brief Set the temperature \f$\mathrm{[K]}\f$ of a fluid phase
     */
    void setTemperature(Scalar value)
    { temperature_ = value; }

protected:
    const FluidState *fs_;
    Scalar temperature_;
};

} // end namespace Dumux

#endif
