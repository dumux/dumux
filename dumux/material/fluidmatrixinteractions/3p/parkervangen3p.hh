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
 * \brief Implementation of van Genuchten's capillary pressure-saturation relation.
 *
 */
#ifndef PARKERVANGEN_3P_HH
#define PARKERVANGEN_3P_HH


#include "parkervangen3pparams.hh"

#include <algorithm>


namespace Dumux
{
/*!
 * \ingroup material
 *
 * \brief Implementation of van Genuchten's capillary pressure <->
 *        saturation relation. This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vince versa.
 *
 * \sa VanGenuchten, VanGenuchtenThreephase
 */
template <class ScalarT, class ParamsT = ParkerVanGen3PParams<ScalarT> >
class ParkerVanGen3P
{

public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief The capillary pressure-saturation curve.
     * \param params Array of parameters
     * \param sw wetting phase saturation
     *
     */
    static Scalar pc(const Params &params, Scalar sw)
    {
        DUNE_THROW(Dune::NotImplemented, "Capillary pressures for three phases is not so simple! Use pcgn, pcnw, and pcgw");
    }
   /*!
     * \brief The capillary pressure-saturation curve copied from MUFTE/pml/constrel3p3cni.c 
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
    this function is just copied from MUFTE/pml/constrel3p3cni.c
    that is why variable names do not yet fulfill Dumux rules, TODO Change */

    Scalar r,se,x,vgm;
    Scalar pc,pcPrime,seRegu;
    Scalar pcvgReg = 0.01;

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
     * \brief The capillary pressure-saturation curve copied from MUFTE/pml/constrel3p3cni.c 
     * \param params Array of parameters
     * \param sw wetting phase saturation or sum of wetting phase saturations
     */
    static Scalar pcnw(const Params &params, Scalar sw)
    {
    /*
         sw = wetting phase saturation, or,
              sum of wetting phase saturations
         alpha : VanGenuchten-alpha
    this function is just copied from MUFTE/pml/constrel3p3cni.c
    that is why variable names do not yet fulfill Dumux rules, TODO Change */

    Scalar r,se,x,vgm;
    Scalar pc,pcPrime,seRegu;
    Scalar pcvgReg = 0.01;

    se   = (sw-params.swr())/(1.-params.snr());

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
                      *(-1/vgm)/params.vgAlpha()/(1-params.snr()-params.swr())/params.vgn();

            /* evaluate tangential */
            r        = (se-seRegu)*pcPrime+pc;
            return(r/params.betaNw());
        }
    }
    /*!
     * \brief The capillary pressure-saturation curve copied from MUFTE/pml/constrel3p3cni.c 
     * \param params Array of parameters
     * \param st sum of wetting (liquid) phase saturations
     */
    static Scalar pcgn(const Params &params, Scalar St)
    {
    /*
         St = sum of wetting (liquid) phase saturations
         alpha : VanGenuchten-alpha
    this function is just copied from MUFTE/pml/constrel3p3cni.c
    that is why variable names do not yet fulfill Dumux rules, TODO Change */

    Scalar r,se,x,vgm;
    Scalar pc,pcPrime,seRegu;
    Scalar pcvgReg = 0.01;

    se   = (St-params.swrx())/(1.-params.swrx());

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
                      *(-1/vgm)/params.vgAlpha()/(1-params.sgr()-params.swrx())/params.vgn();

            /* evaluate tangential */
            r        = (se-seRegu)*pcPrime+pc;
            return(r/params.betaGn());
        }
    }
 /*!
     * \brief The capillary pressure-saturation curve copied from MUFTE/pml/constrel3p3cni.c 
     * \param params Array of parameters
     * \param sn Non-wetting liquid saturation
     */
    static Scalar pcAlpha(const Params &params, Scalar sn)
    {
        /* continuous transition to zero */
        Scalar alpha,sne;

        sne=sn;
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
        DUNE_THROW(Dune::NotImplemented, "sw(pc) for three phases not implemented! Do it yourself!");
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param params Array of parameters
     * \param sw Wetting liquid saturation
    */
    static Scalar dpc_dsw(const Params &params, Scalar sw)
    {
        DUNE_THROW(Dune::NotImplemented, "dpc/dsw for three phases not implemented! Do it yourself!");
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     * \param params Array of parameters
     * \param pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar dsw_dpc(const Params &params, Scalar pc)
    {
        DUNE_THROW(Dune::NotImplemented, "dsw/dpc for three phases not implemented! Do it yourself!");
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
    static Scalar krw(const Params &params,  Scalar saturation, Scalar sn, Scalar sg)
    {

        //transformation to effective saturation
        Scalar se = (saturation - params.swr()) / (1-params.swr());

        /* regularization */
        if(se > 1.0) return 1.;
        if(se < 0.0) return 0.;

        Scalar r = 1. - std::pow(1 - std::pow(se, 1/params.vgm()), params.vgm());
        return std::sqrt(se)*r*r;
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
    static Scalar krn(const Params &params, Scalar sw, Scalar saturation, Scalar sg)
    {

        Scalar swe = std::min((sw - params.swr()) / (1 - params.swr()), 1.);
        Scalar ste = std::min((sw +  saturation - params.swr()) / (1 - params.swr()), 1.);

        // regularization
        if(swe <= 0.0) swe = 0.;
        if(ste <= 0.0) ste = 0.;
        if(ste - swe <= 0.0) return 0.;

        Scalar krn_;
        krn_ = std::pow(1 - std::pow(swe, 1/params.vgm()), params.vgm());
        krn_ -= std::pow(1 - std::pow(ste, 1/params.vgm()), params.vgm());
        krn_ *= krn_;

        if (params.krRegardsSnr())
        {
            // regard Snr in the permeability of the n-phase, see Helmig1997
            Scalar resIncluded = std::max(std::min((saturation - params.snr()/ (1-params.swr())), 1.), 0.);
            krn_ *= std::sqrt(resIncluded );
        }
        else
            krn_ *= std::sqrt(saturation / (1 - params.swr()));   // Hint: (ste - swe) = sn / (1-Srw)


        return krn_;
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
    static Scalar krg(const Params &params, Scalar sw, Scalar sn, Scalar saturation)
    {

        // se = (sw+sn - Sgr)/(1-Sgr)
        Scalar se = std::min(((1-saturation) - params.sgr()) / (1 - params.sgr()), 1.);


        /* regularization */
        if(se > 1.0) return 0.0;
        if(se < 0.0) return 1.0;
        Scalar scalFact = 1.;
        if (saturation<=0.1)
        {
          scalFact = (saturation - params.sgr())/(0.1 - params.sgr());
          if (scalFact < 0.) scalFact = 0.;
        }

        Scalar result = scalFact * std::pow(1 - se, 1.0/3.) * std::pow(1 - std::pow(se, 1/params.vgm()), 2*params.vgm());

        return result;
    }

    /*!
     * \brief The relative permeability for a phase.
     * \param sw Wetting liquid saturation
     * \param sg Gas saturation
     * \param sn Non-wetting liquid saturation
     * \param params Array of parameters.
     * \param phaseIdx indicator, The saturation of all phases.
     */
    static Scalar kr(const Params &params, const int phaseIdx, const Scalar sw, const Scalar sn, const Scalar sg)
    {
        switch (phaseIdx)
        {
        case 0:
            return krw(params, sw, sn, sg);
            break;
        case 1:
            return krn(params, sw, sn, sg);
            break;
        case 2:
            return krg(params, sw, sn, sg);
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
      return params.rhoBulk() * params.KdNAPL();
   }
};
}

#endif
