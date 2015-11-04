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
 * \brief Implementation of a regularized version of van Genuchten's capillary
 *        pressure-saturation relation for three phases.
 */
#ifndef REGULARIZED_PARKERVANGEN_3P_HH
#define REGULARIZED_PARKERVANGEN_3P_HH

#include "parkervangen3p.hh"
#include "regularizedparkervangen3pparams.hh"



#include <dumux/common/spline.hh>

namespace Dumux
{
/*!\ingroup fluidmatrixinteractionslaws
 *
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
    typedef Dumux::ParkerVanGen3P<ScalarT, ParamsT> ParkerVanGen3P;

public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief A regularized Parker- van Genuchten capillary pressure-saturation
     *        curve.
     *
     * regularized part:
     *    - low saturation:  extend the \f$\mathrm{p_c(S_w)}\f$ curve with the slope at the regularization point (i.e. no kink).
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
     * \param sw wetting phase saturation or sum of wetting phase saturations
     *
     */
    static Scalar pcgw(const Params &params, Scalar sw)
    {
    /*
         sw = wetting phase saturation, or,
              sum of wetting phase saturations
         alpha : VanGenuchten-alpha
    */

    Scalar r,se,x,vgm;
    Scalar pc,pcPrime,seRegu;
    Scalar pcvgReg = 0.2;//0.01;

    se   = (sw-params.swr())/(1.-params.sgr());

    /* Snr  = 0.0;   test version   */

    /* regularization */
    if (se<0.0) se=0.0;
    if (se>1.0) se=1.0;
    vgm = 1.-1./params.vgn();

        if (se>pcvgReg && se<1-pcvgReg)
        {
            r = std::pow(se,-1/vgm);
            x = r-1;
            vgm = 1-vgm;
            x = std::pow(x,vgm);
            r = x/params.vgAlpha();
            return(r);
        }
        else
        {
            /* value and derivative at regularization point */
            if (se<=pcvgReg) seRegu = pcvgReg; else seRegu = 1-pcvgReg;
            pc       = std::pow(std::pow(seRegu,-1/vgm)-1,1/params.vgn())/params.vgAlpha();
            pcPrime = std::pow(std::pow(seRegu,-1/vgm)-1,1/params.vgn()-1)*std::pow(seRegu,-1/vgm-1)
                      *(-1/vgm)/params.vgAlpha()/(1-params.sgr()-params.swr())/params.vgn();

            /* evaluate tangential */
            r        = (se-seRegu)*pcPrime+pc;
            return(r/params.betaGw());
        }
    }
  /*!
    
    
    
    
    
    
    
      static Scalar pc(const Params &params, Scalar swe)
    {
        const Scalar sThres = params.thresholdSw();

        // make sure that the capilary pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative is calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (swe <= sThres) {
            Scalar m = BrooksCorey::dpc_dsw(params, sThres);
            Scalar pcsweLow = BrooksCorey::pc(params, sThres);
            return pcsweLow + m*(swe - sThres);
        }
        else if (swe > 1.0) {
            Scalar m = BrooksCorey::dpc_dsw(params, 1.0);
            Scalar pcsweHigh = BrooksCorey::pc(params, 1.0);
            return pcsweHigh + m*(swe - 1.0);
        }

        // if the effective saturation is in an 'reasonable'
        // range, we use the real Brooks-Corey law...
        return BrooksCorey::pc(params, swe);
    }
    
    
    /*!
     * \brief   A regularized Brooks-Corey saturation-capillary pressure curve.
     *
     * regularized part:
     *    - low saturation:  extend the \f$\mathrm{p_c(S_w)}\f$ curve with the slope at the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line (yes, there is a kink :-( ).
     *
     *  The according quantities are obtained by exploiting theorem of intersecting lines.
     *
     * For the non-regularized part:
     *
     * \copydetails BrooksCorey::sw()
     */
    static Scalar sw(const Params &params, Scalar pc)
    {
        const Scalar sThres = params.thresholdSw();

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized version of
        // the Brooks-Corey law
        Scalar swe = BrooksCorey::sw(params, pc);

        // make sure that the capilary pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (swe <= sThres) {
            // invert the low saturation regularization of pc()
            Scalar m = BrooksCorey::dpc_dsw(params, sThres);
            Scalar pcsweLow = BrooksCorey::pc(params, sThres);
            return sThres + (pc - pcsweLow)/m;
        }
        else if (swe > 1.0) {
            Scalar m = BrooksCorey::dpc_dsw(params, 1.0);
            Scalar pcsweHigh = BrooksCorey::pc(params, 1.0);
            return 1.0 + (pc - pcsweHigh)/m;
        }

        return BrooksCorey::sw(params, pc);
    }

    /*!
     * \brief A regularized version of the partial derivative
     *        of the \f$\mathrm{p_c(\overline{S}_w)}\f$ w.r.t. effective saturation
     *        according to Brooks & Corey.
     *
     * regularized part:
     *    - low saturation:  use the slope of the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line and use that slope (yes, there is a kink :-( ).
     *
     * For the non-regularized part:
     *
     * \copydetails BrooksCorey::dpc_dsw()
     */
    static Scalar dpc_dsw(const Params &params, Scalar swe)
    {
        const Scalar sThres = params.thresholdSw();

        // derivative of the regualarization
        if (swe <= sThres) {
            // calculate the slope of the straight line used in pc()
            Scalar m = BrooksCorey::dpc_dsw(params, sThres);
            return m;
        }
        else if (swe > 1.0) {
            // calculate the slope of the straight line used in pc()
            Scalar m = BrooksCorey::dpc_dsw(params, 1.0);
            return m;
        }

        return BrooksCorey::dpc_dsw(params, swe);
    }

    /*!
     * \brief A regularized version of the partial derivative
     *        of the \f$\mathrm{\overline{S}_w(p_c)}\f$ w.r.t. cap.pressure
     *        according to Brooks & Corey.
     *
     *  regularized part:
     *    - low saturation:  use the slope of the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line and use that slope (yes, there is a kink :-( ).
     *
     * For the non-regularized part:
     *
     * \copydetails BrooksCorey::dsw_dpc()
     */
    static Scalar dsw_dpc(const Params &params, Scalar pc)
    {
        const Scalar sThres = params.thresholdSw();

        //instead of return value = inf, return a very large number
        if (params.pe() == 0.0)
        {
            return 1e100;
        }

        // calculate the saturation which corresponds to the
        // saturation in the non-regularized version of the
        // Brooks-Corey law
        Scalar swe;
        if (pc < 0)
            swe = 1.5; // make sure we regularize below
        else
            swe = BrooksCorey::sw(params, pc);

        // derivative of the regularization
        if (swe <= sThres) {
            // calculate the slope of the straight line used in pc()
            Scalar m = BrooksCorey::dpc_dsw(params, sThres);
            return 1/m;
        }
        else if (swe > 1.0) {
            // calculate the slope of the straight line used in pc()
            Scalar m = BrooksCorey::dpc_dsw(params, 1.0);
            return 1/m;
        }
        return 1.0/BrooksCorey::dpc_dsw(params, swe);
    }

    /*!
     * \brief   Regularized version of the  relative permeability
     *          for the wetting phase of
     *          the medium implied by the Brooks-Corey
     *          parameterization.
     *
     *  regularized part:
     *    - below \f$\mathrm{\overline{S}_w =0}\f$:                  set relative permeability to zero
     *    - above \f$\mathrm{\overline{S}_w =1}\f$:                  set relative permeability to one
     *    - between \f$\mathrm{ 0.95 \leq \overline{S}_w \leq 1}\f$:  use a spline as interpolation
     *
     *  For not-regularized part:
        \copydetails BrooksCorey::krw()
     */
    static Scalar krw(const Params &params, Scalar swe)
    {
        if (swe <= 0.0)
            return 0.0;
        else if (swe >= 1.0)
            return 1.0;

        return BrooksCorey::krw(params, swe);
    }

    /*!
     * \brief   Regularized version of the  relative permeability
     *          for the non-wetting phase of
     *          the medium implied by the Brooks-Corey
     *          parameterization.
     *
     * regularized part:
     *    - below \f$\mathrm{\overline{S}_w =0}\f$:                  set relative permeability to zero
     *    - above \f$\mathrm{\overline{S}_w =1}\f$:                  set relative permeability to one
     *    - for \f$\mathrm{0 \leq \overline{S}_w \leq 0.05}\f$:     use a spline as interpolation
     *
         \copydetails BrooksCorey::krn()
     *
     */
    static Scalar krn(const Params &params, Scalar swe)
    {
        if (swe >= 1.0)
            return 0.0;
        else if (swe <= 0.0)
            return 1.0;

        return BrooksCorey::krn(params, swe);
    }
};
}

#endif
