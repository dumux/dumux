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
 *        saturations and takes all other quantities from an other
 *        fluid state.
 */
#ifndef DUMUX_SATURATION_OVERLAY_FLUID_STATE_HH
#define DUMUX_SATURATION_OVERLAY_FLUID_STATE_HH

namespace Dumux {

/*!
 * \ingroup FluidStates
 * \brief This is a fluid state which allows to set the fluid
 *        saturations and takes all other quantities from an other
 *        fluid state.
 */
template <class FluidState>
class SaturationOverlayFluidState
{
public:
    static constexpr int numPhases = FluidState::numPhases;
    static constexpr int numComponents = FluidState::numComponents;

    //! export the scalar type
    using Scalar =  typename FluidState::Scalar;

    /*!
     * \brief Constructor
     *
     * \param fs Fluidstate
     * The overlay fluid state copies the saturation from the argument,
     * so it initially behaves exactly like the underlying fluid
     * state.
     */
    SaturationOverlayFluidState(const FluidState &fs)
    : fs_(&fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            saturation_[phaseIdx] = fs.saturation(phaseIdx);
    }

    // copy & move constructor / assignment operators
    SaturationOverlayFluidState(const SaturationOverlayFluidState &fs) = default;
    SaturationOverlayFluidState(SaturationOverlayFluidState &&fs) = default;
    SaturationOverlayFluidState& operator=(const SaturationOverlayFluidState &fs) = default;
    SaturationOverlayFluidState& operator=(SaturationOverlayFluidState &&fs) = default;

    /*****************************************************
     * Generic access to fluid properties (No assumptions
     * on thermodynamic equilibrium required)
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
     *  @copydoc CompositionalFluidState::moleFraction()
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    { return fs_->moleFraction(phaseIdx, compIdx); }

    /*!
     *  @copydoc CompositionalFluidState::massFraction()
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
     *  @copydoc CompositionalFluidState::temperature()
     */
    Scalar temperature(int phaseIdx) const
    { return fs_->temperature(phaseIdx); }

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
     * \brief Set the saturation \f$\mathrm{[-]}\f$ of a fluid phase
     */
    void setSaturation(int phaseIdx, Scalar value)
    { saturation_[phaseIdx] = value; }

protected:
    const FluidState *fs_;
    Scalar saturation_[numPhases] = {};
};

} // end namespace Dumux

#endif
