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
 * \brief   Parameters that are necessary for the \em regularization of
 *          the Parker - Van Genuchten capillary pressure model.
 */
#ifndef DUMUX_REGULARIZED_PARKERVANGEN_3P_PARAMS_HH
#define DUMUX_REGULARIZED_PARKERVANGEN_3P_PARAMS_HH

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

#include "parkervangen3pparams.hh"

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief   Parameters that are necessary for the \em regularization of
 *          the Parker - van Genuchten capillary pressure model.
 */
template <class ScalarT>
class RegularizedParkerVanGen3PParams : public ParkerVanGen3PParams<ScalarT>
{
    using ParkerVanGen3PParams = Dumux::ParkerVanGen3PParams<ScalarT>;

public:
    using Scalar = ScalarT;

    RegularizedParkerVanGen3PParams()
        : ParkerVanGen3PParams(), constRegularization_(false)
    {
        pcLowS_ = 1e-2;
        pcHighS_ = 99e-2;
    }

    RegularizedParkerVanGen3PParams(Scalar vgAlpha, Scalar vgn, Scalar KdNAPL, Scalar rhoBulk,
                         Dune::FieldVector<Scalar, 4> residualSaturation, Scalar betaNw = 1.,
                         Scalar betaGn = 1., Scalar betaGw = 1., bool regardSnr=false)
        : ParkerVanGen3PParams(vgAlpha, vgn, KdNAPL, rhoBulk,
                         residualSaturation, betaNw,
                         betaGn , betaGw, regardSnr), constRegularization_(false)
    {
        pcLowS_ = 1e-2;
        pcHighS_ = 99e-2;
    }

    /*!
     * \brief Threshold saturation below which the capillary pressure
     *        is regularized.
     *
     * This is just 1%. If you need a different value, overload this
     * class.
     */
    Scalar pcLowS() const
    {
        // Most problems are very sensitive to this value
        // (e.g. making it smaller might result in negative
        // pressures)
        //
        // If you want to use a different regularization threshold,
        // overload this class and supply the new class as second
        // template parameter for the RegularizedVanGenuchten law!
        return pcLowS_;
    }

    /*!
     * \brief Threshold saturation above which the capillary pressure
     *        is regularized.
     *
     * This is just 99%. If you need a different value, overload this
     * class.
     */
    Scalar pcHighS() const
    {
        // Most problems are very sensitive to this value
        // (e.g. making it smaller might result in negative
        // pressures)
        //
        // If you want to use a different regularization threshold,
        // overload this class and supply the new class as second
        // template parameter for the RegularizedVanGenuchten law!
        return pcHighS_;
    }

    /*!
     * \brief Set the lower saturation threshold value
     * \param input The saturation threshold value
     */
    void setPcLowS(const Scalar input)
    { pcLowS_ = input; }

    /*!
     * \brief Set the upper saturation threshold value
     * \param input The saturation threshold value
     */
    void setPcHighS(const Scalar input)
    { pcHighS_ = input; }

    /*!
     * \brief Choose whether to use a constant value for regularization of the
     *        pc-S curves or not
     * \param input True or false
     */
    void useConstRegularization(const bool input)
    {
        constRegularization_ = input;
    }

    /*!
     * \brief Returns whether to use a constant value for regularization of the
     *        pc-S curves or not
     */
    bool constRegularization() const
    {
        return constRegularization_;
    }

private:
    Scalar pcLowS_;
    Scalar pcHighS_;
    bool constRegularization_;

};
} // namespace Dumux

#endif
