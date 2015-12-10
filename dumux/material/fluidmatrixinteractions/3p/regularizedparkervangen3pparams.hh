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
 * \brief   Parameters that are necessary for the \em regularization of
 *          the Parker - Van Genuchten capillary pressure model.
 */
#ifndef DUMUX_REGULARIZED_PARKERVANGEN_3P_PARAMS_HH
#define DUMUX_REGULARIZED_PARKERVANGEN_3P_PARAMS_HH

#include "parkervangen3pparams.hh"

namespace Dumux
{
/*!
 * \brief   Parameters that are necessary for the \em regularization of
 *          the Parker - van Genuchten capillary pressure model.
 *
 *        \ingroup fluidmatrixinteractionsparams
 */
template <class ScalarT>
class RegularizedParkerVanGen3PParams : public Dumux::ParkerVanGen3PParams<ScalarT>
{
    typedef Dumux::ParkerVanGen3PParams<ScalarT> ParkerVanGen3PParams;

public:
    typedef ScalarT Scalar;

    RegularizedParkerVanGen3PParams()
        : ParkerVanGen3PParams()
    { thresholdSw_ = 1e-2; }

    RegularizedParkerVanGen3PParams(Scalar vgAlpha, Scalar vgn, Scalar KdNAPL, Scalar rhoBulk,
                         Dune::FieldVector<Scalar, 4> residualSaturation, Scalar betaNw = 1.,
                         Scalar betaGn = 1., Scalar betaGw = 1., bool regardSnr=false)
        : ParkerVanGen3PParams(vgAlpha, vgn, KdNAPL, rhoBulk,
                         residualSaturation, betaNw,
                         betaGn , betaGw, regardSnr)
    { thresholdSw_ = 1e-2; }

    /*!
     * \brief Threshold saturation below which the capillary pressure
     *        is regularized.
     *
     * This is just 1%. If you need a different value, overload this
     * class.
     */
    Scalar thresholdSw() const
    {
        // Most problems are very sensitive to this value
        // (e.g. making it smaller might result in negative
        // pressures)
        //
        // If you want to use a different regularization threshold,
        // overload this class and supply the new class as second
        // template parameter for the regularized Parker - van Genuchten law!
        return thresholdSw_;
    }

    /*!
     * \brief Set the saturation threshold value
     * \param input The saturation threshold value
     */
    void setThresholdSw(const Scalar input)
    { thresholdSw_ = input; }

private:
    Scalar thresholdSw_;

};
} // namespace Dumux

#endif
