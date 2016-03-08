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
 * \brief This is a fluid state which allows to set the fluid
 *        pressures and takes all other quantities from an other
 *        fluid state.
 */
#ifndef DUMUX_PRESSURE_OVERLAY_FLUID_STATE_HH
#define DUMUX_PRESSURE_OVERLAY_FLUID_STATE_HH

#include <array>

#include <dumux/common/valgrind.hh>

namespace Dumux
{
/*!
 * \ingroup FluidStates
 * \brief This is a fluid state which allows to set the fluid
 *        pressures and takes all other quantities from an other
 *        fluid state.
 */
template <class Scalar, class FluidState>
class PressureOverlayFluidState
{
public:
    enum { numPhases = FluidState::numPhases };
    enum { numComponents = FluidState::numComponents };

    /*!
     * \brief Constructor
     *
     * \param fs Fluidstate
     * The overlay fluid state copies the pressures from the argument,
     * so it initially behaves exactly like the underlying fluid
     * state.
     */
    PressureOverlayFluidState(const FluidState &fs)
        : fs_(&fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            pressure_[phaseIdx] = fs.pressure(phaseIdx);
    }

    // copy constructor
    PressureOverlayFluidState(const PressureOverlayFluidState &fs)
        : fs_(fs.fs_)
        , pressure_(fs.pressure_)
    {
    }

    // assignment operator
    PressureOverlayFluidState &operator=(const PressureOverlayFluidState &fs)
    {
        fs_ = fs.fs_;
        pressure_ = fs.pressure_;
        return *this;
    }

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
     * @copydoc CompositionalFluidState::moleFraction()
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    { return fs_->moleFraction(phaseIdx, compIdx); }

    /*!
     * @copydoc CompositionalFluidState::massFraction()
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    { return fs_->massFraction(phaseIdx, compIdx); }

    /*!
     * @copydoc CompositionalFluidState::averageMolarMass()
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return fs_->averageMolarMass(phaseIdx); }

    /*!
     * @copydoc CompositionalFluidState::molarity()
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    { return fs_->molarity(phaseIdx, compIdx); }

    /*!
     * @copydoc CompositionalFluidState::fugacity()
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    { return fs_->fugacity(phaseIdx, compIdx); }

    /*!
     * @copydoc CompositionalFluidState::fugacityCoefficient()
     */
    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    { return fs_->fugacityCoefficient(phaseIdx, compIdx); }

    /*!
     * @copydoc CompositionalFluidState::molarVolume()
     */
    Scalar molarVolume(int phaseIdx) const
    { return fs_->molarVolume(phaseIdx); }

    /*!
     * @copydoc CompositionalFluidState::density()
     */
    Scalar density(int phaseIdx) const
    { return fs_->density(phaseIdx); }

    /*!
     *@copydoc CompositionalFluidState::molarDensity()
     */
    Scalar molarDensity(int phaseIdx) const
    { return fs_->molarDensity(phaseIdx); }

    /*!
     * @copydoc CompositionalFluidState::temperature()
     */
    Scalar temperature(int phaseIdx) const
    { return fs_->temperature(phaseIdx); }

    /*!
     *  @copydoc CompositionalFluidState::pressure()
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; }

    /*!
     * @copydoc CompositionalFluidState::enthalpy()
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
     * \brief Set the pressure \f$\mathrm{[Pa]}\f$ of a fluid phase
     */
    void setPressure(int phaseIdx, Scalar value)
    { pressure_[phaseIdx] = value; }

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
        Valgrind::CheckDefined(pressure_);
    }

protected:
    const FluidState *fs_;
    std::array<Scalar, numPhases> pressure_;
};

} // end namespace Dumux

#endif
