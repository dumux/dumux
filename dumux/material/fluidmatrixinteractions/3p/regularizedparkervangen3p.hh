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
 * \brief Implementation of a regularized version of van Genuchten's capillary
 *        pressure-saturation relation for three phases.
 */
#ifndef REGULARIZED_PARKERVANGEN_3P_HH
#define REGULARIZED_PARKERVANGEN_3P_HH

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

#include <algorithm>

#include "parkervangen3p.hh"
#include "regularizedparkervangen3pparams.hh"

#include <dumux/common/spline.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the regularized van Genuchten's
 *        capillary pressure <-> saturation relation.
 *        This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vince versa.
 *
 *        In order to avoid very steep gradients the marginal values are "regularized".
 *        This means that in stead of following the curve of the material law in these regions, some linear approximation is used.
 *        Doing this is not worse than following the material law. E.g. for very low wetting phase values the material
 *        laws predict infinite values for \f$\mathrm{p_c}\f$ which is completely unphysical. In case of very high wetting phase
 *        saturations the difference between regularized and "pure" material law is not big.
 *
 *        Regularizing has the additional benefit of being numerically friendly: Newton's method does not like infinite gradients.
 *
 *        The implementation is accomplished as follows:
 *        - check whether we are in the range of regularization
 *         - yes: use the regularization
 *         - no: forward to the standard material law.
 *
 *         For an example figure of the regularization: RegularizedVanGenuchten
 *
 * \see BrooksCorey
 */
template <class ScalarT, class ParamsT = RegularizedParkerVanGen3PParams<ScalarT> >
class RegularizedParkerVanGen3P
{
    using ParkerVanGen3P = Dumux::ParkerVanGen3P<ScalarT, ParamsT>;

public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    /*!
     * \brief A regularized Parker- van Genuchten capillary pressure-saturation
     *        curve.
     *
     * regularized part:
     *    - low saturation: extend the \f$\mathrm{p_c(S_w)}\f$ curve with the slope at the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line (yes, there is a kink :-( ).
     *
     * For the non-regularized part:
     *
     * \copydetails ParkerVanGen3P::pc()
     */
    static Scalar pc(const Params &params, Scalar sw)
    {
        return ParkerVanGen3P::pc(params, sw);
    }

    /*!
     * \brief The capillary pressure-saturation curve for the gas and wetting phase
     * \param params Array of parameters
     * \param swe Effective wetting phase saturation
     *
     */
    static Scalar pcgw(const Params &params, Scalar swe)
    {
        //if specified, a constant value is used for regularization
        using std::clamp;
        if (params.constRegularization())
            swe = clamp(swe, 0.0, 1.0);

        // retrieve the low and the high threshold saturations for the
        // unregularized capillary pressure curve from the parameters
        const Scalar swThLow = params.pcLowS();
        const Scalar swThHigh = params.pcHighS();

        // make sure that the capillary pressure observes a derivative
        // != 0 for 'illegal' saturations. This is favourable for the
        // newton solver (if the derivative is calculated numerically)
        // in order to get the saturation moving to the right
        // direction if it temporarily is in an 'illegal' range.
        if (swe < swThLow)
        {
            const Scalar mLow = ParkerVanGen3P::dpcgw_dswe(params, swThLow);
            return ParkerVanGen3P::pcgw(params, swThLow) + mLow*(swe - swThLow);
        }
        else if (swe > swThHigh)
        {
            Scalar yTh = ParkerVanGen3P::pcgw(params, swThHigh);
            Scalar m1 = (0.0 - yTh)/(1.0 - swThHigh)*2;

            if (swe < 1.0) {
                // use spline between threshold swe and 1.0
                Scalar mTh = ParkerVanGen3P::dpcgw_dswe(params, swThHigh);
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
        return ParkerVanGen3P::pcgw(params, swe);
    }

    /*!
     * \brief The capillary pressure-saturation curve for the non-wettigng and wetting phase
     * \param params Array of parameters
     * \param swe Effective wetting phase saturation
     */
    static Scalar pcnw(const Params &params, Scalar swe)
    {
        //if specified, a constant value is used for regularization
        if(params.constRegularization())
        {
            if(swe < 0.0)
                swe = 0.0;
            if(swe > 1.0)
                swe = 1.0;
        }

        // retrieve the low and the high threshold saturations for the
        // unregularized capillary pressure curve from the parameters
        const Scalar swThLow = params.pcLowS();
        const Scalar swThHigh = params.pcHighS();

        // make sure that the capillary pressure observes a derivative
        // != 0 for 'illegal' saturations. This is favourable for the
        // newton solver (if the derivative is calculated numerically)
        // in order to get the saturation moving to the right
        // direction if it temporarily is in an 'illegal' range.
        if (swe < swThLow)
        {
            const Scalar mLow = ParkerVanGen3P::dpcnw_dswe(params, swThLow);
            return ParkerVanGen3P::pcnw(params, swThLow) + mLow*(swe - swThLow);
        }
        else if (swe > swThHigh)
        {
            Scalar yTh = ParkerVanGen3P::pcnw(params, swThHigh);
            Scalar m1 = (0.0 - yTh)/(1.0 - swThHigh)*2;

            if (swe < 1.0) {
                // use spline between threshold swe and 1.0
                Scalar mTh = ParkerVanGen3P::dpcnw_dswe(params, swThHigh);
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
        return ParkerVanGen3P::pcnw(params, swe);
    }

    /*!
     * \brief The capillary pressure-saturation curve for the gas and nonwetting phase
     * \param params Array of parameters
     * \param ste Effective total liquid (wetting + nonwetting) saturation
     */
    static Scalar pcgn(const Params &params, Scalar ste)
    {
        //if specified, a constant value is used for regularization
        if(params.constRegularization())
        {
            if(ste < 0.0)
                ste = 0.0;
            if(ste > 1.0)
                ste = 1.0;
        }

        // retrieve the low and the high threshold saturations for the
        // unregularized capillary pressure curve from the parameters
        const Scalar swThLow = params.pcLowS();
        const Scalar swThHigh = params.pcHighS();

        // make sure that the capillary pressure observes a derivative
        // != 0 for 'illegal' saturations. This is favourable for the
        // newton solver (if the derivative is calculated numerically)
        // in order to get the saturation moving to the right
        // direction if it temporarily is in an 'illegal' range.
        if (ste < swThLow)
        {
            const Scalar mLow = ParkerVanGen3P::dpcgn_dste(params, swThLow);
            return ParkerVanGen3P::pcgn(params, swThLow) + mLow*(ste - swThLow);
        }
        else if (ste > swThHigh)
        {
            Scalar yTh = ParkerVanGen3P::pcgn(params, swThHigh);
            Scalar m1 = (0.0 - yTh)/(1.0 - swThHigh)*2;

            if (ste < 1.0) {
                // use spline between threshold swe and 1.0
                Scalar mTh = ParkerVanGen3P::dpcgn_dste(params, swThHigh);
                Spline<Scalar> sp(swThHigh, 1.0, // x0, x1
                                  yTh, 0, // y0, y1
                                  mTh, m1); // m0, m1
                return sp.eval(ste);
            }
            else {
                // straight line for swe > 1.0
                return m1*(ste - 1.0) + 0.0;
            }
        }

        // if the effective saturation is in an 'reasonable'
        // range, we use the real van genuchten law...
        return ParkerVanGen3P::pcgn(params, ste);
    }

    /*!
     * \brief This function ensures a continuous transition from 2 to 3 phases and vice versa
     * \param params Array of parameters
     * \param sne Effective nonwetting liquid saturation
     */
    static Scalar pcAlpha(const Params &params, Scalar sne)
    {
        return ParkerVanGen3P::pcAlpha(params, sne);
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     * \param params Array of parameters
     * \param pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar sw(const Params &params, Scalar pc)
    {
        return ParkerVanGen3P::sw(params, pc);
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param params Array of parameters
     * \param swe Effective wetting liquid saturation
     */
    static Scalar dpc_dswe(const Params &params, Scalar swe)
    {
        return ParkerVanGen3P::dpc_dswe(params, swe);
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     * \param params Array of parameters
     * \param pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar dswe_dpc(const Params &params, Scalar pc)
    {
        return ParkerVanGen3P::dswe_dpc(params, pc);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by van Genuchten's
     *        parameterization.
     *
     * The permeability of water in a 3p system equals the standard 2p description.
     * (see p61. in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83.)
     *
     * \param params Array of parameters.
     * \param swe Effective wetting phase saturation
     */
    static Scalar krw(const Params &params,  const Scalar swe)
    {
        //use regularization
        if(swe > 1.0) return 1.0;
        if(swe < 0.0) return 0.0;

        //or use actual material law
        return ParkerVanGen3P::krw(params, swe);
    }

    /*!
     * \brief The relative permeability for the nonwetting phase
     *        after the Model of Parker et al. (1987).
     *
     * See model 7 in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83.
     * or more comprehensive in
     * "Estimation of primary drainage three-phase relative permeability for organic
     * liquid transport in the vadose zone", Leonardo I. Oliveira, Avery H. Demond,
     * Journal of Contaminant Hydrology 66 (2003), 261-285
     *
     * \param params Array of parameters.
     * \param swe Effective wetting phase saturation
     * \param sn Absolute nonwetting liquid saturation
     * \param ste Effective total liquid (wetting + nonwetting) saturation
     */
    static Scalar krn(const Params &params, Scalar swe, Scalar sn, Scalar ste)
    {
        using std::clamp;
        swe = clamp(swe, 0.0, 1.0);
        ste = clamp(ste, 0.0, 1.0);

        if(ste - swe <= 0.0)
            return 0.0;

        //or use actual material law
        return ParkerVanGen3P::krn(params, swe, sn, ste);
    }


    /*!
     * \brief The relative permeability for the nonwetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * The permeability of gas in a 3p system equals the standard 2p description.
     * (see p61. in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83.)
     *
     * \param params Array of parameters.
     * \param ste Effective total liquid (wetting + nonwetting) saturation
     */
    static Scalar krg(const Params &params, const Scalar ste)
    {
        //return 0 if there is no gas
        if(ste > 1.0)
            return 0.0;

        // use linear regularization for very high gas saturations
        // to avoid a kink in the curve and to maintain a slope for
        // the Newton solver
        const Scalar threshold = 1e-3;
        if(ste <= threshold)
        {
            const Scalar mSwr = ParkerVanGen3P::dkrg_dste(params, threshold);
            const Scalar ySwr = ParkerVanGen3P::krg(params, threshold);
            return ySwr + mSwr*(ste - threshold);
        }

        // For very low gas saturations:
        // We use a scaling factor that decreases the gas phase permeability quite fast at very low gas phase
        // saturations, thus making that phase virtually immobile.
        // This prevents numerical issues related to the degeneration of the gas phase mass balance for the 3p3c model
        // at very low gas phase saturations.

        //get the absolute gas phase saturation
        const Scalar st = ste*(1 - params.swr()) + params.swr();
        const Scalar sg = 1.0 - st;
        using std::max;
        const Scalar scalFact = (sg > 0.1) ? 1.0 : max(0.0,
                                                      (sg - params.sgr())/(0.1 - params.sgr()));

        //or use actual material law
        return scalFact * ParkerVanGen3P::krg(params, ste);
    }

    /*!
     * \brief The relative permeability for a phase.
     * \param params Array of parameters.
     * \param phaseIdx Indicator, The saturation of all phases.
     * \param swe Effective wetting phase saturation
     * \param sn Absolute nonwetting liquid saturation
     * \param ste Effective total liquid (wetting + nonwetting) saturation
     */
    static Scalar kr(const Params &params, const int phaseIdx, const Scalar swe, const Scalar sn, const Scalar ste)
    {
        switch (phaseIdx)
        {
        case 0:
            return krw(params, swe);
        case 1:
            return krn(params, swe, sn, ste);
        case 2:
            return krg(params, ste);
        }
        DUNE_THROW(Dune::InvalidStateException,
                   "Invalid phase index ");
    }

   /*!
    * \brief the basis for calculating adsorbed NAPL in storage term
    * \param params Array of parameters
    */
   static Scalar bulkDensTimesAdsorpCoeff (const Params &params)
   {
      return ParkerVanGen3P::bulkDensTimesAdsorpCoeff(params);
   }

};
} // end namespace Dumux

#endif
