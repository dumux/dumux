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
 * \brief Specification of the material params for the van Genuchten
 *        capillary pressure model.
 *
 * In comparison to the 2p version, this parameter container also includes
 * the residual saturations, as their inclusion is very model-specific.
 */
#ifndef PARKERVANGEN_PARAMS_3P_HH
#define PARKERVANGEN_PARAMS_3P_HH

#include <dune/common/fvector.hh>

namespace Dumux
{
/*!
 * \brief Reference implementation of a van Genuchten params
 */
template<class ScalarT>
class ParkerVanGen3PParams
{
public:
    typedef ScalarT Scalar;

    ParkerVanGen3PParams()
    {betaGw_ = betaNw_ = betaGn_ = 1.;}

    ParkerVanGen3PParams(Scalar vgAlpha, Scalar vgn, Scalar KdNAPL, Scalar rhoBulk,
                         Dune::FieldVector<Scalar, 4> residualSaturation, Scalar betaNw = 1.,
                         Scalar betaGn = 1., Scalar betaGw = 1., bool regardSnr=false)
    {
        setVgAlpha(vgAlpha);
        setVgn(vgn);
        setSwr(residualSaturation[0]);
        setSnr(residualSaturation[1]);
        setSgr(residualSaturation[2]);
        setSwrx(residualSaturation[3]);
        setKrRegardsSnr(regardSnr);
        setKdNAPL(KdNAPL);
        setBetaNw(betaNw);
        setBetaGn(betaGn);
        setBetaGw(betaGw);
        setRhoBulk(rhoBulk);
    }

    /*!
     * \brief Return the \f$\alpha\f$ shape parameter of van Genuchten's
     *        curve.
     */
    Scalar vgAlpha() const
    { return vgAlpha_; }

    /*!
     * \brief Set the \f$\alpha\f$ shape parameter of van Genuchten's
     *        curve.
     */
    void setVgAlpha(Scalar v)
    { vgAlpha_ = v; }

    /*!
     * \brief Return the \f$m\f$ shape parameter of van Genuchten's
     *        curve.
     */
    Scalar vgm() const
    { return vgm_; }

    /*!
     * \brief Set the \f$m\f$ shape parameter of van Genuchten's
     *        curve.
     *
     * The \f$n\f$ shape parameter is set to \f$n = \frac{1}{1 - m}\f$
     */
    void setVgm(Scalar m)
    { vgm_ = m; vgn_ = 1/(1 - vgm_); }

    /*!
     * \brief Return the \f$n\f$ shape parameter of van Genuchten's
     *        curve.
     */
    Scalar vgn() const
    { return vgn_; }

    /*!
     * \brief Set the \f$n\f$ shape parameter of van Genuchten's
     *        curve.
     *
     * The \f$n\f$ shape parameter is set to \f$m = 1 - \frac{1}{n}\f$
     */
    void setVgn(Scalar n)
    { vgn_ = n; vgm_ = 1 - 1/vgn_; }

    /*!
     * \brief Return the residual saturation.
     */
    Scalar satResidual(int phaseIdx) const
    {
        switch (phaseIdx)
        {
        case 0:
            return swr_;
            break;
        case 1:
            return snr_;
            break;
        case 2:
            return sgr_;
            break;
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Set all residual saturations.
     */
    void setResiduals(Dune::FieldVector<Scalar, 3> residualSaturation)
    {
        setSwr(residualSaturation[0]);
        setSnr(residualSaturation[1]);
        setSgr(residualSaturation[2]);
    }


    /*!
     * \brief Return the residual wetting saturation.
     */
    Scalar swr() const
    { return swr_; }

    /*!
     * \brief Set the residual wetting saturation.
     */
    void setSwr(Scalar input)
    { swr_ = input; }

    /*!
     * \brief Return the residual non-wetting saturation.
     */
    Scalar snr() const
    { return snr_; }

    /*!
     * \brief Set the residual non-wetting saturation.
     */
    void setSnr(Scalar input)
    { snr_ = input; }

    /*!
     * \brief Return the residual gas saturation.
     */
    Scalar sgr() const
    { return sgr_; }

    /*!
     * \brief Set the residual gas saturation.
     */
    void setSgr(Scalar input)
    { sgr_ = input; }

    Scalar swrx() const
    { return swrx_; }

    /*!
     * \brief Set the residual gas saturation.
     */
    void setSwrx(Scalar input)
    { swrx_ = input; }

    /*!
     * \brief defines the scaling parameters of capillary pressure between the phases (=1 for Gas-Water)
     */
    void setBetaNw(Scalar input)
    { betaNw_ = input; }

    void setBetaGn(Scalar input)
    { betaGn_ = input; }

    void setBetaGw(Scalar input)
    { betaGw_ = input; }

    /*!
     * \brief Return the values for the beta scaling parameters of capillary pressure between the phases
     */
    Scalar betaNw() const
    { return betaNw_; }

    Scalar betaGn() const
    { return betaGn_; }

    Scalar betaGw() const
    { return betaGw_; }

    /*!
     * \brief defines if residual n-phase saturation should be regarded in its relative permeability.
     */
    void setKrRegardsSnr(bool input)
    { krRegardsSnr_ = input; }

    /*!
     * \brief Calls if residual n-phase saturation should be regarded in its relative permeability.
     */
    bool krRegardsSnr() const
    { return krRegardsSnr_; }


    /*!
     * \brief Return the bulk density of the porous medium
     */
    Scalar rhoBulk() const
    { return rhoBulk_; }

    /*!
     * \brief Set the bulk density of the porous medium
     */
    void setRhoBulk(Scalar input)
    { rhoBulk_ = input; }

    /*!
     * \brief Return the adsorption coefficient
     */
    Scalar KdNAPL() const
    { return KdNAPL_; }

    /*!
     * \brief Set the adsorption coefficient
     */
    void setKdNAPL(Scalar input)
    { KdNAPL_ = input; }


private:
    Scalar vgAlpha_;
    Scalar vgm_;
    Scalar vgn_;
    Scalar swr_;
    Scalar snr_;
    Scalar sgr_;
    Scalar swrx_;     /* (sw+sn)_r */

    Scalar KdNAPL_;
    Scalar rhoBulk_;

    Scalar betaNw_;
    Scalar betaGn_;
    Scalar betaGw_;

    bool krRegardsSnr_ ;
};
} // namespace Dumux

#endif
