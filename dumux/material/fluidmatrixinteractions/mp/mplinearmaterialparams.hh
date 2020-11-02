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
 * \ingroup Fluidmatrixinteractions
 * \brief Reference implementation of parameters for the M-phase linear
 * material material.
 */
#ifndef MP_LINEAR_MATERIAL_PARAMS_HH
#define MP_LINEAR_MATERIAL_PARAMS_HH

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Reference implementation of params for the linear M-phase
 *        material material.
 */
template<int numPhasesV, class ScalarT>
class MpLinearMaterialParams
{
public:
    using Scalar = ScalarT;
    enum { numPhases = numPhasesV };

    [[deprecated("Use FluidMatrix::MPLinearMaterial. Class will be removed after 3.3.")]]
    MpLinearMaterialParams()
    {
        for (int i = 0; i < numPhases; ++i) {
            setPcMinSat(i, 0.0);
            setPcMaxSat(i, 0.0);
        }
    }

    /*!
     * \brief Return the threshold saturation at which the relative
     *        permeability starts to get regularized.
     * \param phaseIdx Index of the phase
     *
     * This is simply 10%
     */
    Scalar sReg(int phaseIdx) const
    { return 0.10; }

    /*!
     * \brief Return the capillary pressure in \f$\mathrm{[Pa]}\f$ for a phase \f$\mathrm{\alpha}\f$ at \f$\mathrm{S_\alpha=0}\f$.
     * \param phaseIdx Index of the phase
     */
    Scalar pcMinSat(int phaseIdx) const
    { return pcMinSat_[phaseIdx]; }

    /*!
     * \brief Set the capillary pressure in \f$\mathrm{[Pa]}\f$ for a phase \f$\mathrm{\alpha}\f$ at \f$\mathrm{S_\alpha=0}\f$.
     * \param phaseIdx Index of the phase
     * \param val Value of the capillary pressure
     */
    void setPcMinSat(int phaseIdx, Scalar val)
    { pcMinSat_[phaseIdx] = val; }

    /*!
     * \brief Return the capillary pressure in \f$\mathrm{[Pa]}\f$ for a phase \f$\mathrm{\alpha}\f$ at \f$\mathrm{S_\alpha=1}\f$.
     * \param phaseIdx Index of the phase
     */
    Scalar pcMaxSat(int phaseIdx) const
    { return pcMaxSat_[phaseIdx]; }

    /*!
     * \brief Set the capillary pressure in \f$\mathrm{[Pa]}\f$ for a phase \f$\mathrm{\alpha}\f$ at \f$\mathrm{S_\alpha=1}\f$.
     * \param phaseIdx Index of the phase
     * \param val Value of the capillary pressure
     */
    void setPcMaxSat(int phaseIdx, Scalar val)
    { pcMaxSat_[phaseIdx] = val; }

    /*!
     * \brief Return the threshold saturation respective phase below
     *        which the relative permeability gets regularized.
     *
     * This is just 5%. If you need a different value, write your own
     * parameter class.
     * \param phaseIdx Index of the phase
     */
    Scalar krLowS(int phaseIdx) const
    { return 0.05; }

    /*!
     * \brief Return the threshold saturation of the respective phase
     *        above which the relative permeability gets regularized.
     *
     * This is just 95%. If you need a different value, write your own
     * parameter class.
     * \param phaseIdx Index of the phase
     */
    Scalar krHighS(int phaseIdx) const
    { return 0.95; }

private:
    Scalar pcMaxSat_[numPhases];
    Scalar pcMinSat_[numPhases];
};
} // namespace Dumux

#endif
