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
    static Scalar pcgw(const Params &params, Scalar swe)
    {
        Scalar r,x,vgm;
        Scalar pc,pcPrime,seRegu;
        Scalar pcvgReg = 0.2;//0.01; //TODO paramter


        /* regularization */
        if (swe<0.0) swe=0.0;
        if (swe>1.0) swe=1.0;
        vgm = 1.-1./params.vgn();

        if (swe>pcvgReg && swe<1-pcvgReg) //use actual material law
        {
            ParkerVanGen3P::pcgw(params, swe);
        }
        else //use regularization
        {
            /* value and derivative at regularization point */
            if (swe<=pcvgReg) seRegu = pcvgReg; else seRegu = 1-pcvgReg;
            pc       = std::pow(std::pow(seRegu,-1/vgm)-1,1/params.vgn())/params.vgAlpha();
            pcPrime = std::pow(std::pow(seRegu,-1/vgm)-1,1/params.vgn()-1)*std::pow(seRegu,-1/vgm-1)
                      *(-1/vgm)/params.vgAlpha()/(1-params.sgr()-params.swr())/params.vgn();

            /* evaluate tangential */
            r        = (swe-seRegu)*pcPrime+pc;
            return(r/params.betaGw());
        }
    }

  /*!
     * \brief The capillary pressure-saturation curve for the non-wettigng and wetting phase
     * \param params Array of parameters
     * \param sw wetting phase saturation or sum of wetting phase saturations
     */
    static Scalar pcnw(const Params &params, Scalar swe)
    {

    Scalar r,vgm;
    Scalar pc,pcPrime,seRegu;
    Scalar pcvgReg = /*0.2*/0.01;

    /* regularization */
    if (swe<0.0) swe=0.0;
    if (swe>1.0) swe=1.0;
    vgm = 1.-1./params.vgn();

        if (swe>pcvgReg && swe<1-pcvgReg)
        {
          return ParkerVanGen3P::pcnw(params, swe);
        }
        else
        {//TODO: snr??
            /* value and derivative at regularization point */
            if (swe<=pcvgReg) seRegu = pcvgReg; else seRegu = 1-pcvgReg;
            pc       = std::pow(std::pow(seRegu,-1/vgm)-1,1/params.vgn())/params.vgAlpha();
            pcPrime = std::pow(std::pow(seRegu,-1/vgm)-1,1/params.vgn()-1)*std::pow(seRegu,-1/vgm-1)
                      *(-1/vgm)/params.vgAlpha()/(1-params.snr()-params.swr())/params.vgn();

            /* evaluate tangential */
            r        = (swe-seRegu)*pcPrime+pc;
            return(r/params.betaNw());
        }
    }
    /*!
     * \brief The capillary pressure-saturation curve for the gas and non-wetting phase
     * \param params Array of parameters
     * \param St sum of wetting (liquid) phase saturations
     */
    static Scalar pcgn(const Params &params, Scalar ste)
    {
    Scalar r,vgm;
    Scalar pc,pcPrime,seRegu;
    Scalar pcvgReg = /*0.2*/0.01;


    /* regularization */
    if (ste<0.0) ste=0.0;
    if (ste>1.0) ste=1.0;
    vgm = 1.-1./params.vgn();

        if (ste>pcvgReg && ste<1-pcvgReg)
        {
          return ParkerVanGen3P::pcgn(params, ste);
        }
        else
        {//TODO: swrx
            /* value and derivative at regularization point */
            if (ste<=pcvgReg) seRegu = pcvgReg; else seRegu = 1-pcvgReg;
            pc       = std::pow(std::pow(seRegu,-1/vgm)-1,1/params.vgn())/params.vgAlpha();
            pcPrime = std::pow(std::pow(seRegu,-1/vgm)-1,1/params.vgn()-1)*std::pow(seRegu,-1/vgm-1)
                      *(-1/vgm)/params.vgAlpha()/(1-params.sgr()-params.swrx())/params.vgn();

            /* evaluate tangential */
            r        = (ste-seRegu)*pcPrime+pc;
            return(r/params.betaGn());
        }
    }
    /*!
     * \brief The capillary pressure-saturation curve copied from MUFTE/pml/constrel3p3cni.c
     * \param params Array of parameters
     * \param sn Non-wetting liquid saturation
     */
    static Scalar pcAlpha(const Params &params, Scalar sne)
    {
        /* continuous transition to zero */
        Scalar alpha;
        //TODO: gehoert das hier hin??

        /* regularization */
        if (sne<=0.001) sne=0.0;
        if (sne>=1.0) sne=1.0;

        if (sne>params.snr()) alpha = 1.0;
        else
        {
         if (params.snr()>=0.001) alpha = sne/params.snr();
         else          alpha = 0.0;
        }
        return(alpha);
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
     * \param sw Wetting liquid saturation
    */
    static Scalar dpc_dsw(const Params &params, Scalar sw)
    {
        return ParkerVanGen3P::dpc_dsw(params, sw);
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     * \param params Array of parameters
     * \param pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar dsw_dpc(const Params &params, Scalar pc)
    {
        return ParkerVanGen3P::dsw_dpc(params, pc);
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
     * \param sn Non-wetting liquid saturation
     * \param sg Gas saturation
     * \param saturation wetting liquid saturation
     * \param params Array of parameters.
     */
    static Scalar krw(const Params &params,  Scalar swe)
    {
        /* regularization */
        if(swe > 1.0) return 1.;
        if(swe < 0.0) return 0.;

        return ParkerVanGen3P::krw(params, swe);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        after the Model of Parker et al. (1987).
     *
     * See model 7 in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83.
     * or more comprehensive in
     * "Estimation of primary drainage three-phase relative permeability for organic
     * liquid transport in the vadose zone", Leonardo I. Oliveira, Avery H. Demond,
     * Journal of Contaminant Hydrology 66 (2003), 261-285
     *
     *
     * \param sw Wetting liquid saturation
     * \param sg Gas saturation
     * \param saturation Non-wetting liquid saturation
     * \param params Array of parameters.
     */
    static Scalar krn(const Params &params, Scalar swe, Scalar sne, Scalar ste)
    {
        swe = std::min(swe, 1.);
        ste = std::min(ste, 1.);

        // regularization
        if(swe <= 0.0) swe = 0.;
        if(ste <= 0.0) ste = 0.;
        if(ste - swe <= 0.0) return 0.;

        return ParkerVanGen3P::krn(params, swe, sne, ste);
    }


    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * The permeability of gas in a 3p system equals the standard 2p description.
     * (see p61. in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83.)
     *
     * \param sw Wetting liquid saturation
     * \param sn Non-wetting liquid saturation
     * \param saturation Gas saturation
     * \param params Array of parameters.
     */
    static Scalar krg(const Params &params, Scalar ste)
    {

        /* regularization */
        if(ste > 1.0) return 0.0;
        if(ste < 0.0) return 1.0;

        return ParkerVanGen3P::krg(params, ste);
    }

    /*!
     * \brief The relative permeability for a phase.
     * \param sw Wetting liquid saturation
     * \param sg Gas saturation
     * \param sn Non-wetting liquid saturation
     * \param params Array of parameters.
     * \param phaseIdx indicator, The saturation of all phases.
     */ //TODO: indices???
    static Scalar kr(const Params &params, const int phaseIdx, const Scalar swe, const Scalar sne, const Scalar ste/*sg*/)
    {
        switch (phaseIdx)
        {
        case 0:
            return krw(params, swe);
            break;
        case 1:
            return krn(params, swe, sne, ste);
            break;
        case 2:
            return krg(params, ste);
            break;
        }
        return 0;
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
}

#endif
