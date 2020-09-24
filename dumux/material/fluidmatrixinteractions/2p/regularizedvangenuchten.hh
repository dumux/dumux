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
 * \brief Implementation of the regularized version of the van Genuchten's
 *        capillary pressure / relative permeability  <-> saturation relation.
 */
#ifndef REGULARIZED_VAN_GENUCHTEN_HH
#define REGULARIZED_VAN_GENUCHTEN_HH

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

#include "vangenuchten.hh"
#include "regularizedvangenuchtenparams.hh"

#include <algorithm>

#include <dumux/common/spline.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the regularized  van Genuchten's
 *        capillary pressure / relative permeability  <-> saturation relation.
 *
 *        This class bundles the "raw" curves as
 *        static members and doesn't concern itself converting
 *        absolute to effective saturations and vice versa.
 *
 *        In order to avoid very steep gradients the marginal values
 *        are "regularized".  This means that in stead of following
 *        the curve of the material law in these regions, some linear
 *        approximation is used.  Doing this is not worse than
 *        following the material law. E.g. for very low wetting phase
 *        values the material laws predict infinite values for
 *        \f$\mathrm{p_c}\f$ which is completely unphysical. In case of very
 *        high wetting phase saturations the difference between
 *        regularized and "pure" material law is not big.
 *
 *        Regularizing has the additional benefit of being numerically
 *        friendly: Newton's method does not like infinite gradients.
 *
 *        The implementation is accomplished as follows:
 *        - check whether we are in the range of regularization
 *         - yes: use the regularization
 *         - no: forward to the standard material law.
 *
 *        An example of the regularization of the capillary pressure curve is shown below:
 *        \image html regularizedVanGenuchten.png
 *
 * \see VanGenuchten
 */
template <class ScalarT, class ParamsT = RegularizedVanGenuchtenParams<ScalarT> >
class RegularizedVanGenuchten
{
    using VanGenuchten = Dumux::VanGenuchten<ScalarT, ParamsT>;

public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    /*!
     * \brief A regularized van Genuchten capillary pressure-saturation
     *          curve.
     *
     * regularized part:
     *    - low saturation:  extend the \f$\mathrm{p_c(S_w)}\f$ curve with the slope at the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line (yes, there is a kink :-( ).
     *
     *  For not-regularized part:
     *
     *   \copydetails VanGenuchten::pc()
     */
    static Scalar pc(const Params &params, Scalar swe)
    {
        // retrieve the low and the high threshold saturations for the
        // unregularized capillary pressure curve from the parameters
        const Scalar swThLow = params.pcLowSw();
        const Scalar swThHigh = params.pcHighSw();

        // make sure that the capillary pressure observes a derivative
        // != 0 for 'illegal' saturations. This is favourable for the
        // newton solver (if the derivative is calculated numerically)
        // in order to get the saturation moving to the right
        // direction if it temporarily is in an 'illegal' range.
        if (swe < swThLow) {
            return VanGenuchten::pc(params, swThLow) + mLow_(params)*(swe - swThLow);
        }
        else if (swe > swThHigh)
        {
            Scalar yTh = VanGenuchten::pc(params, swThHigh);
            Scalar m1 = (0.0 - yTh)/(1.0 - swThHigh)*2;

            if (swe < 1.0) {
                // use spline between threshold swe and 1.0
                Scalar mTh = VanGenuchten::dpc_dswe(params, swThHigh);
                Spline<Scalar> sp(swThHigh, 1.0, // x0, x1
                                  yTh, 0, // y0, y1
                                  mTh, m1); // m0, m1
                return sp.eval(swe);
            }
            else {
                // straight line for swe > 1.0
                return m1*(swe - 1.0) + 0.0;
            }
        }

        // if the effective saturation is in an 'reasonable'
        // range, we use the real van genuchten law...
        return VanGenuchten::pc(params, swe);
    }

    /*!
     * \brief   A regularized van Genuchten saturation-capillary pressure curve.
     *
     * regularized part:
     *    - low saturation:  extend the \f$\mathrm{p_c(S_w)}\f$ curve with the slope at the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line (yes, there is a kink :-( ).
     *
     *  The according quantities are obtained by exploiting theorem of intersecting lines.
     *
     *  For not-regularized part:
     *
     *    \copydetails VanGenuchten::sw()
     *
     */
    static Scalar sw(const Params &params, Scalar pc)
    {
        // retrieve the low and the high threshold saturations for the
        // unregularized capillary pressure curve from the parameters
        const Scalar swThLow = params.pcLowSw();
        const Scalar swThHigh = params.pcHighSw();

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized verision of van
        // Genuchten's law
        Scalar sw;
        if (pc <= 0) {
            // for swThHigh = 1.0 the slope would get infinity
            // swThHigh > 1.0 are not sensible threshold values
            // setting swThHigh = 1.0 is a way to disable regularization
            if (swThHigh > 1.0 - std::numeric_limits<Scalar>::epsilon())
                return 1.0;

            // invert straight line for swe > 1.0
            Scalar yTh = VanGenuchten::pc(params, swThHigh);
            Scalar m1 = (0.0 - yTh)/(1.0 - swThHigh)*2;
            return pc/m1 + 1.0;
        }
        else
            sw = VanGenuchten::sw(params, pc);

        // invert the regularization if necessary
        if (sw <= swThLow) {
            // invert the low saturation regularization of pc()
            Scalar pcswLow = VanGenuchten::pc(params, swThLow);
            return (pc - pcswLow)/mLow_(params) + swThLow;
        }
        else if (sw > swThHigh)
        {
            Scalar yTh = VanGenuchten::pc(params, swThHigh);
            Scalar m1 = (0.0 - yTh)/(1.0 - swThHigh)*2;

            // invert spline between threshold swe and 1.0
            Scalar mTh = VanGenuchten::dpc_dswe(params, swThHigh);
            Spline<Scalar> sp(swThHigh, 1.0, // x0, x1
                              yTh, 0, // m0, m1
                              mTh, m1); // m0, m1
            return sp.intersectInterval(swThHigh, 1.0,
                                        0, 0, 0, pc);
        }

        return sw;
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar endPointPc(const Params &params)
    { return 0.0; }

    /*!
     * \brief A regularized version of the partial derivative
     *        of the \f$\mathrm{p_c(\overline{S}_w)}\f$ w.r.t. effective saturation
     *        according to van Genuchten.
     *
     * regularized part:
     *    - low saturation:  use the slope of the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line and use that slope (yes, there is a kink :-( ).
     *
     *        For not-regularized part:
     *
     * \copydetails VanGenuchten::dpc_dswe()
     */
    static Scalar dpc_dswe(const Params &params, Scalar swe)
    {
        // retrieve the low and the high threshold saturations for the
        // unregularized capillary pressure curve from the parameters
        const Scalar swThLow = params.pcLowSw();
        const Scalar swThHigh = params.pcHighSw();

        // derivative of the regularization
        if (swe < swThLow) {
            // the slope of the straight line used in pc()
            return mLow_(params);
        }
        else if (swe > swThHigh)
        {
            Scalar yTh = VanGenuchten::pc(params, swThHigh);
            Scalar m1 = (0.0 - yTh)/(1.0 - swThHigh)*2;

            if (swe < 1.0) {
                // use spline between threshold swe and 1.0
                Scalar mTh = VanGenuchten::dpc_dswe(params, swThHigh);
                Spline<Scalar> sp(swThHigh, 1.0, // x0, x1
                                  yTh, 0, // y0, y1
                                  mTh, m1); // m0, m1
                return sp.evalDerivative(swe);
            }
            else {
                // straight line for swe > 1.0
                return m1;
            }
        }

        return VanGenuchten::dpc_dswe(params, swe);
    }

    /*!
     * \brief A regularized version of the partial derivative
     *        of the \f$\mathrm{\overline{S}_w(p_c)}\f$ w.r.t. cap.pressure
     *        according to van Genuchten.
     *
     *  regularized part:
     *    - low saturation:  use the slope of the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line and use that slope (yes, there is a kink :-( ).
     *
     *        For not-regularized part:
     * \copydetails VanGenuchten::dswe_dpc()
     */
    static Scalar dswe_dpc(const Params &params, Scalar pc)
    {
        const Scalar swThLow = params.pcLowSw();
        const Scalar swThHigh = params.pcHighSw();

        if (pc <= 0)
        {
            // for swThHigh = 1.0 the slope gets infinity
            // swThHigh > 1.0 are not sensible threshold values
            // setting swThHigh = 1.0 is a way to disable regularization
            // return the maximum representable number with the given precision
            if (swThHigh > 1.0 - std::numeric_limits<Scalar>::epsilon())
                return std::numeric_limits<Scalar>::max();

            Scalar yTh = VanGenuchten::pc(params, swThHigh);
            Scalar m1 = (0.0 - yTh)/(1.0 - swThHigh)*2;
            return 1.0/m1;
        }

        const auto sw = VanGenuchten::sw(params, pc);

        // derivative of the regularization
        if (sw <= swThLow)
            return 1/mLow_(params);

        if (sw > swThHigh)
        {
            const Scalar yTh = VanGenuchten::pc(params, swThHigh);
            const Scalar m1 = (0.0 - yTh)/(1.0 - swThHigh)*2;

            // invert spline between threshold swe and 1.0
            const Scalar mTh = VanGenuchten::dpc_dswe(params, swThHigh);
            const Spline<Scalar> sp(swThHigh, 1.0, // x0, x1
                                    yTh, 0, // m0, m1
                                    mTh, m1); // m0, m1
            const auto swReg = sp.intersectInterval(swThHigh, 1.0,
                                                    0, 0, 0, pc);
            // derivative of the inverse of the function is one over derivative of the function
            return 1.0/sp.evalDerivative(swReg);
        }

        return VanGenuchten::dswe_dpc(params, pc);
    }

    /*!
     * \brief   Regularized version of the  relative permeability
     *          for the wetting phase of
     *          the medium implied by the van Genuchten
     *          parameterization.
     *
     *  regularized part:
     *    - below \f$\mathrm{\overline{S}_w =0}\f$:                  set relative permeability to zero
     *    - above \f$\mathrm{\overline{S}_w =1}\f$:                  set relative permeability to one
     *    - between \f$\mathrm{0.95 \leq \overline{S}_w \leq 1}\f$:  use a spline as interpolation
     *
     *  For not-regularized part:
     * \copydetails VanGenuchten::krw()
     */
    static Scalar krw(const Params &params, Scalar swe)
    {
        // retrieve the high threshold saturation for the
        // unregularized relative permeability curve of the wetting
        // phase from the parameters
        const Scalar swThHigh = params.krwHighSw();

        if (swe < 0)
            return 0;
        else if (swe > 1 - std::numeric_limits<Scalar>::epsilon())
            return 1;
        else if (swe > swThHigh) {
            using Spline = Dumux::Spline<Scalar>;
            Spline sp(swThHigh, 1.0, // x1, x2
                      VanGenuchten::krw(params, swThHigh), 1.0, // y1, y2
                      VanGenuchten::dkrw_dswe(params, swThHigh), 0); // m1, m2
            return sp.eval(swe);
        }

        return VanGenuchten::krw(params, swe);
    }

    /*!
     * \brief A regularized version of the derivative of the relative
     *        permeability for the wetting phase in regard to the wetting
     *        saturation of the medium implied by the van Genuchten parameterization.
     *
     * \copydetails VanGenuchten::dkrw_dswe()
     */
    static Scalar dkrw_dswe(const Params &params, Scalar swe)
    {
        // retrieve the high threshold saturation for the
        // unregularized relative permeability curve of the wetting
        // phase from the parameters
        const Scalar swThHigh = params.krwHighSw();

        // derivative of the regularization
        // the slope is zero below sw=0.0 and above sw=1.0
        if (swe < 0)
            return 0.0;
        else if (swe > 1 - std::numeric_limits<Scalar>::epsilon())
            return 0.0;

        // for high sw we need the slope of the interpolation spline
        else if (swe > swThHigh) {
            using Spline = Dumux::Spline<Scalar>;
            Spline sp(swThHigh, 1.0, // x1, x2
                      VanGenuchten::krw(params, swThHigh), 1.0, // y1, y2
                      VanGenuchten::dkrw_dswe(params, swThHigh), 0); // m1, m2
            return sp.evalDerivative(swe);
        }

        return VanGenuchten::dkrw_dswe(params, swe);
    }

    /*!
     * \brief   Regularized version of the  relative permeability
     *          for the nonwetting phase of
     *          the medium implied by the van Genuchten
     *          parameterization.
     *
     * regularized part:
     *    - below \f$\mathrm{\overline{S}_w =0}\f$:                  set relative permeability to zero
     *    - above \f$\mathrm{\overline{S}_w =1}\f$:                  set relative permeability to one
     *    - for \f$\mathrm{0 \leq \overline{S}_w \leq 0.05}\f$:     use a spline as interpolation
     *
     * \copydetails VanGenuchten::krn()
     *
     */
    static Scalar krn(const Params &params, Scalar swe)
    {
        // retrieve the low threshold saturation for the unregularized
        // relative permeability curve of the nonwetting phase from
        // the parameters
        const Scalar swThLow = params.krnLowSw();

        if (swe <= 0)
            return 1;
        else if (swe >= 1)
            return 0;
        else if (swe < swThLow) {
            using Spline = Dumux::Spline<Scalar>;
            Spline sp(0.0, swThLow, // x1, x2
                      1.0, VanGenuchten::krn(params, swThLow), // y1, y2
                      0.0, VanGenuchten::dkrn_dswe(params, swThLow)); // m1, m2
            return sp.eval(swe);
        }

        return VanGenuchten::krn(params, swe);
    }

    /*!
     * \brief A regularized version of the derivative of the relative permeability
     *        for the nonwetting phase in regard to the wetting saturation of
     *        the medium as implied by the van Genuchten parameterization.
     *
     * \copydetails VanGenuchten::dkrw_dswe()
     */
    static Scalar dkrn_dswe(const Params &params, Scalar swe)
    {
        // retrieve the low threshold saturation for the unregularized
        // relative permeability curve of the nonwetting phase from
        // the parameters
        const Scalar swThLow = params.krnLowSw();

        if (swe <= 0)
            return 0.0;
        else if (swe >= 1)
            return 0.0;
        else if (swe < swThLow) {
            using Spline = Dumux::Spline<Scalar>;
            Spline sp(0.0, swThLow, // x1, x2
                      1.0, VanGenuchten::krn(params, swThLow), // y1, y2
                      0.0, VanGenuchten::dkrn_dswe(params, swThLow)); // m1, m2
            return sp.evalDerivative(swe);
        }

        return VanGenuchten::dkrn_dswe(params, swe);
    }

private:
    // the slope of the straight line used to regularize saturations
    // below the minimum saturation

    /*!
     * \brief   The slope of the straight line used to regularize
     *          saturations below the minimum saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar mLow_(const Params &params)
    {
        const Scalar swThLow = params.pcLowSw();

        return VanGenuchten::dpc_dswe(params, swThLow);
    }

    /*!
     * \brief   The slope of the straight line used to regularize
     *          saturations above the minimum saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar mHigh_(const Params &params)
    {
        const Scalar swThHigh = params.pcHighSw();

        // for swThHigh = 1.0 the slope would get infinity
        // swThHigh > 1.0 are not sensible threshold values
        // setting swThHigh = 1.0 is a way to disable regularization
        if (swThHigh > 1.0 - std::numeric_limits<Scalar>::epsilon())
            return 0.0;

        Scalar pcswHigh = VanGenuchten::pc(params, swThHigh);
        return (0 - pcswHigh)/(1.0 - swThHigh);
    }
};

} // end namespace Dumux

#endif
