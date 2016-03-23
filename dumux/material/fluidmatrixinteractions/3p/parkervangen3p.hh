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
 * \brief Implementation of van Genuchten's capillary pressure-saturation relation for three phases.
 *
 */
#ifndef PARKERVANGEN_3P_HH
#define PARKERVANGEN_3P_HH


#include "parkervangen3pparams.hh"

#include <algorithm>

#warning The Parker-VanGenchten 3P material law \
         has been thoroughly revised. The conversion from \
         absolute to effective saturations and regularization \
         are now done in additional separate classes: \
         <dumux/material/fluidmatrixinteractions/3p/efftoabslaw.hh> and \
         <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3p.hh> \
         Make sure to use these classes in your spatialParams. \
         This warning will be removed after the next release of DuMux.


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
    static Scalar pc(const Params &params, const Scalar sw)
    {
        DUNE_THROW(Dune::NotImplemented, "Capillary pressures for three phases is not so simple! Use pcgn, pcnw, and pcgw");
    }
   /*!
     * \brief The capillary pressure-saturation curve for the gas and wetting phase
     * \param params Array of parameters
     * \param swe Effective wetting phase saturation
     *
     */
    static Scalar pcgw(const Params &params, const Scalar swe)
    {
        assert(0 <= swe && swe <= 1);
        return pc_(params, swe);
    }

  /*!
     * \brief The capillary pressure-saturation curve for the non-wettigng and wetting phase
     * \param params Array of parameters
     * \param swe Effective wetting phase saturation
     */
    static Scalar pcnw(const Params &params, const Scalar swe)
    {
        assert(0 <= swe && swe <= 1);
        return pc_(params, swe)/params.betaNw();
    }

    /*!
     * \brief The capillary pressure-saturation curve for the gas and non-wetting phase
     * \param params Array of parameters
     * \param ste Effective total liquid (wetting + non-wetting) saturation
     */
    static Scalar pcgn(const Params &params, const Scalar ste)
    {
        assert(0 <= ste && ste <= 1);
        return pc_(params,ste)/params.betaGn();
    }

     /*!
     * \brief This function ensures a continous transition from 2 to 3 phases and vice versa
     * \param params Array of parameters
     * \param sne Non-wetting liquid saturation
     */
    static Scalar pcAlpha(const Params &params, Scalar sne)
    {
        Scalar alpha;

        /* regularization */
        if (sne<=0.001) sne=0.0;
        if (sne>=1.0) sne=1.0;

        if (sne>params.snr())
            alpha = 1.0;
        else
        {
            if (params.snr()>=0.001)
                alpha = sne/params.snr();
            else
                alpha = 0.0;
        }
        return alpha;
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     * \param params Array of parameters
     * \param pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar sw(const Params &params, const Scalar pc)
    {
        DUNE_THROW(Dune::NotImplemented, "sw(pc) for three phases not implemented! Do it yourself!");
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param params Array of parameters
     * \param swe Effective wetting liquid saturation
    */
    DUNE_DEPRECATED_MSG("dpc_dsw(const Params &params, const Scalar swe) is deprecated. Use dpc_dswe(const Params &params, const Scalar swe) instead.")
    static Scalar dpc_dsw(const Params &params, const Scalar swe)
    {
        return dpc_dswe(params, swe);
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param params Array of parameters
     * \param swe Effective wetting liquid saturation
    */
    static Scalar dpc_dswe(const Params &params, const Scalar swe)
    {
        DUNE_THROW(Dune::NotImplemented, "dpc/dswe for three phases not implemented! Do it yourself!");
    }

     /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param params Array of parameters
     * \param seRegu Effective wetting phase saturation for regularization
    */
    DUNE_DEPRECATED_MSG("dpcgw_dsw(const Params &params, const Scalar seRegu) is deprecated. Use dpcgw_dswe(const Params &params, const Scalar seRegu) instead.")
    static Scalar dpcgw_dsw(const Params &params, const Scalar seRegu)
    {
        return dpcgw_dswe(params, seRegu);
    }

     /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param params Array of parameters
     * \param seRegu Effective wetting phase saturation for regularization
    */
    static Scalar dpcgw_dswe(const Params &params, const Scalar seRegu)
    {
        const Scalar powSeRegu = pow(seRegu, -1/params.vgm());
        return - 1.0/params.vgAlpha() * pow(powSeRegu - 1, 1.0/params.vgn() - 1)/params.vgn()
            * powSeRegu/seRegu/params.vgm()/params.betaGw();
    }

     /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param params Array of parameters
     * \param seRegu Effective wetting phase saturation for regularization
    */
    DUNE_DEPRECATED_MSG("dpcnw_dsw(const Params &params, const Scalar seRegu) is deprecated. Use dpcnw_dswe(const Params &params, const Scalar seRegu) instead.")
    static Scalar dpcnw_dsw(const Params &params, const Scalar seRegu)
    {
        return dpcnw_dswe(params, seRegu);
    }

     /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param params Array of parameters
     * \param seRegu Effective wetting phase saturation for regularization
    */
    static Scalar dpcnw_dswe(const Params &params, const Scalar seRegu)
    {
        const Scalar powSeRegu = pow(seRegu, -1/params.vgm());
        return - 1.0/params.vgAlpha() * pow(powSeRegu - 1, 1.0/params.vgn() - 1)/params.vgn()
            * powSeRegu/seRegu/params.vgm()/params.betaNw();
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param params Array of parameters
     * \param seRegu Effective wetting phase saturation for regularization
    */
    DUNE_DEPRECATED_MSG("dpcgn_dst(const Params &params, const Scalar seRegu) is deprecated. Use dpcgn_dste(const Params &params, const Scalar seRegu) instead.")
    static Scalar dpcgn_dst(const Params &params, const Scalar seRegu)
    {
        return dpcgn_dste(params, seRegu);
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param params Array of parameters
     * \param seRegu Effective wetting phase saturation for regularization
    */
    static Scalar dpcgn_dste(const Params &params, const Scalar seRegu)
    {
        const Scalar powSeRegu = pow(seRegu, -1/params.vgm());
        return - 1.0/params.vgAlpha() * pow(powSeRegu - 1, 1.0/params.vgn() - 1)/params.vgn()
            * powSeRegu/seRegu/params.vgm()/params.betaGn();
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     * \param params Array of parameters
     * \param pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    DUNE_DEPRECATED_MSG("dsw_dpc(const Params &params, const Scalar pc) is deprecated. Use dswe_dpc(const Params &params, const Scalar pc) instead.")
    static Scalar dsw_dpc(const Params &params, const Scalar pc)
    {
        return dswe_dpc(params, pc);
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     * \param params Array of parameters
     * \param pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar dswe_dpc(const Params &params, const Scalar pc)
    {
        DUNE_THROW(Dune::NotImplemented, "dswe/dpc for three phases not implemented! Do it yourself!");
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by van Genuchten's
     *        parameterization.
     *
     * The permeability of water in a 3p system equals the standard 2p description.
     * (see p61. in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83.) \cite delshad1989 <BR>
     *
     * \param params Array of parameters.
     * \param swe Effective wetting phase saturation
     */
    static Scalar krw(const Params &params,  const Scalar swe)
    {
        const Scalar r = 1.0 - std::pow(1 - std::pow(swe, 1/params.vgm()), params.vgm());
        return std::sqrt(swe)*r*r;
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        after the Model of Parker et al. (1987).
     *
     * See model 7 in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83 \cite delshad1989 <BR>
     * or more comprehensive in
     * "Estimation of primary drainage three-phase relative permeability for organic
     * liquid transport in the vadose zone", Leonardo I. Oliveira, Avery H. Demond,
     * Journal of Contaminant Hydrology 66 (2003), 261-285 \cite oliveira2003 <BR>
     *
     *
     * \param params Array of parameters.
     * \param swe Effective wetting phase saturation
     * \param sn Absolute non-wetting liquid saturation
     * \param ste Effective total liquid (wetting + non-wetting) saturation
     */
    static Scalar krn(const Params &params, const Scalar swe, const Scalar sn, const Scalar ste)
    {
        Scalar krn;
        krn = std::pow(1 - std::pow(swe, 1/params.vgm()), params.vgm());
        krn -= std::pow(1 - std::pow(ste, 1/params.vgm()), params.vgm());
        krn *= krn;

        if (params.krRegardsSnr())
        {
            // regard Snr in the permeability of the n-phase, see Helmig1997
            Scalar resIncluded = std::max(std::min((sn - params.snr()/ (1-params.swr())), 1.0), 0.0);
            krn *= std::sqrt(resIncluded );
        }
        else
            krn *= std::sqrt(sn / (1 - params.swr()));   // Hint: (ste - swe) = sn / (1-Srw)

        return krn;
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * The permeability of gas in a 3p system equals the standard 2p description.
     * (see p61. in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83.) \cite delshad1989 <BR>
     *
     * \param params Array of parameters.
     * \param ste Effective total liquid (wetting + non-wetting) saturation
     */
    static Scalar krg(const Params &params, const Scalar ste)
    {
        assert(0 <= ste && ste <= 1);
        return std::cbrt(1 - ste) * std::pow(1 - std::pow(ste, 1/params.vgm()), 2*params.vgm());
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        gas phase in regard to the total liquid saturation of
     *        the medium as implied by the van Genuchten
     *        parameterization.
     *
     * \param ste The mobile total liquid saturation.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar dkrg_dste(const Params &params, Scalar ste)
    {
        assert(0 < ste && ste <= 1);

        const Scalar x = std::pow(ste, 1.0/params.vgm());
        return
            -std::pow(1.0 - x, 2*params.vgm())
            *std::pow(1.0 - ste, -2.0/3)
            *(1.0/3 + 2*x/ste);
    }

    /*!
     * \brief The relative permeability for a phase.
     * \param params Array of parameters.
     * \param phaseIdx indicator, The saturation of all phases.
     * \param swe Effective wetting phase saturation
     * \param sn Absolute non-wetting liquid saturation
     * \param ste Effective total liquid (wetting + non-wetting) saturation
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
      return params.rhoBulk() * params.KdNAPL();
   }

private:

    /*!
     * \brief The standard van Genuchten two-phase pc-S relation either with respect to
     *        the effective wetting phase saturation Swe or the effective total liquid saturation Ste.
     * \param params Array of parameters.
     * \param Se Effective wetting phase ortotal liquid  saturation
     */
    const static Scalar pc_(const Params &params, const Scalar se)
    {
        return std::pow(std::pow(se, -1/params.vgm()) - 1, 1/params.vgn())/params.vgAlpha();
    }

};
}

#endif
