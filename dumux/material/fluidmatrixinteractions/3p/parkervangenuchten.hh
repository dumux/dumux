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
 * \brief Implementation of van Genuchten's capillary pressure-saturation relation for three phases.
 */
#ifndef PARKER_VANGENUCHTEN_3P_HH
#define PARKER_VANGENUCHTEN_3P_HH

#include <algorithm>
#include <dumux/common/optionalscalar.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/spline.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/material/fluidmatrixinteractions/2p/noregularization.hh>

namespace Dumux::FluidMatrix {

struct ParkerVanGenuchten3PEffToAbsPolicy
{
    /*!
     * \brief The parameter type
     * \tparam Scalar The scalar type
     * \note The efftoabs policy need two parameters: \f$\mathrm{S_{w,r}}, \mathrm{S_{n,r}}\f$.
     *       For the respective formulas check out the description of the free function.
     */
    template<class Scalar>
    struct Params
    {
        Params(const Scalar swr = 0.0, const Scalar snr = 0.0, const Scalar sgr = 0.0)
        : swr_(swr), snr_(snr), sgr_(sgr)
        {}

        /*!
         * \brief Return the residual wetting saturation.
         */
        Scalar swr() const
        { return swr_; }

        /*!
         * \brief Set the residual wetting saturation.
         */
        void setSwr(Scalar v)
        { swr_ = v; }

        /*!
         * \brief Return the residual nonwetting saturation.
         */
        Scalar snr() const
        { return snr_; }

        /*!
         * \brief Set the residual nonwetting saturation.
         */
        void setSnr(Scalar v)
        { snr_ = v; }
        /*!

         * \brief Return the residual gas phase saturation.
         */
        Scalar sgr() const
        { return sgr_; }

        /*!
         * \brief Set the residual gas phase saturation.
         */
        void setSgr(Scalar v)
        { sgr_ = v; }

        bool operator== (const Params& p) const
        {
            return Dune::FloatCmp::eq(swr(), p.swr(), 1e-6)
                   && Dune::FloatCmp::eq(snr(), p.snr(), 1e-6)
                   && Dune::FloatCmp::eq(sgr(), p.sgr(), 1e-6);
        }
    private:
        Scalar swr_;
        Scalar snr_;
        Scalar sgr_;
    };

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    template<class Scalar>
    static Params<Scalar> makeParams(const std::string& paramGroup)
    {
        Params<Scalar> params;
        params.setSwr(getParamFromGroup<Scalar>(paramGroup, "Swr", 0.0));
        params.setSnr(getParamFromGroup<Scalar>(paramGroup, "Snr", 0.0));
        params.setSgr(getParamFromGroup<Scalar>(paramGroup, "Sgr", 0.0));
        return params;
    }

    /*!
     * \brief Convert an absolute wetting saturation to an effective one.
     *
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{[{S}_w]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective saturation of the wetting phase.
     */
    template<class Scalar>
    static Scalar swToSwe(const Scalar sw, const Params<Scalar>& params)
    {
        return (sw - params.swr())/(1.0 - params.swr()); // TODO other residual saturations?
    }

    /*!
     * \brief Convert an effective wetting saturation to an absolute one.
     *
     * \param swe Effective saturation of the nonwetting phase \f$\mathrm{[\overline{S}_n]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Absolute saturation of the nonwetting phase.
     */
    template<class Scalar>
    static Scalar sweToSw(const Scalar swe, const Params<Scalar>& params)
    {
        return swe*(1.0 - params.swr()) + params.swr();  // TODO other residual saturations?
    }

    /*!
     * \brief Derivative of the effective saturation w.r.t. the absolute saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the effective saturation w.r.t. the absolute saturation.
     */
    template<class Scalar>
    static Scalar dswe_dsw(const Params<Scalar>& params)
    {
        return 1.0/(1.0 - params.swr()); // TODO other residual saturations?
    }

    /*!
     * \brief Derivative of the absolute saturation w.r.t. the effective saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the absolute saturation w.r.t. the effective saturation.
     */
    template<class Scalar>
    static Scalar dsw_dswe(const Params<Scalar>& params)
    {
        return 1.0 - params.swr(); // TODO other residual saturations?
    }

    /*!
     * \brief Convert an absolute nonwetting saturation to an effective one.
     *
     * \param sn Absolute saturation of the nonwetting phase \f$\mathrm{[{S}_n]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective saturation of the nonwetting phase.
     */
    template<class Scalar>
    static Scalar snToSne(const Scalar sn, const Params<Scalar>& params)
    {
        return sn; // sne equals sn  // TODO other residual saturations?
    }

    /*!
     * \brief Convert an absolute total liquid saturation to an effective one.
     *
     * \param st Absolute saturation of the total liquid phase (sw+sn) \f$\mathrm{[{S}_n]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective saturation of the nonwetting phase.
     */
    template<class Scalar>
    static Scalar stToSte(const Scalar st, const Params<Scalar>& params)
    {
        return (st-params.swr()) / (1-params.swr()); // TODO other residual saturations?
    }

    /*!
     * \brief Convert an effective wetting saturation to an absolute one.
     *
     * \param ste Effective total liquid (wetting + nonwetting) saturation
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Absolute saturation of the nonwetting phase.
     */
    template<class Scalar>
    static Scalar steToSt(const Scalar ste, const Params<Scalar>& params)
    {
        return ste*(1.0 - params.swr()) + params.swr();  // TODO other residual saturations?
    }

    /*!
     * \brief Derivative of the effective saturation w.r.t. the absolute saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the effective saturation w.r.t. the absolute saturation.
     */
    template<class Scalar>
    static Scalar dste_dst(const Params<Scalar>& params)
    {
        return 1.0/(1.0 - params.swr() /*- params.snr() - params.sgr()*/); // TODO other residual saturations?
    }

    /*!
     * \brief Derivative of the absolute saturation w.r.t. the effective saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the absolute saturation w.r.t. the effective saturation.
     */
    template<class Scalar>
    static Scalar dst_dste(const Params<Scalar>& params)
    {
        return 1.0 - params.swr(); // TODO other residual saturations?
    }
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of Parker/vanGenuchten's capillary pressure <->
 *        saturation relation for three phases. This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vince versa.
 */
class ParkerVanGenuchten3P
{

public:
    /*!
     * \brief The parameter type
     * \tparam Scalar The scalar type
     * \note The Parker/vanGenuchten laws are parameterized with four parameters: \f$\mathrm{n, m, \alpha, l}\f$.
     *
     * - \f$\mathrm{\alpha}\f$ shape parameter \f$\mathrm{[1/Pa]}\f$
     * - \f$\mathrm{n}\f$ shape parameter \f$\mathrm{[-]}\f$
     * - \f$\mathrm{swr}\f$ wetting phase residual saturation \f$\mathrm{[-]}\f$
     * - \f$\mathrm{swr}\f$ nonwetting phase residual saturation \f$\mathrm{[-]}\f$
     * - \f$\mathrm{betaNw}\f$ scaling parameter \f$\mathrm{[-]}\f$
     * - \f$\mathrm{betaGn}\f$ scaling parameter \f$\mathrm{[-]}\f$
     * - \f$\mathrm{betaGw}\f$ scaling parameter \f$\mathrm{[-]}\f$
     * - \f$\mathrm{regardSnr}\f$ determines whether snr is considered for krn or not
     */
    template<class Scalar>
    struct Params
    {
        Params(Scalar alpha, Scalar n, Scalar swr = 0.0, Scalar snr = 0.0,
               Scalar betaNw = 1.0, Scalar betaGn = 1.0, Scalar betaGw = 1.0, bool regardSnr = false)
        : alpha_(alpha), n_(n), m_(1.0 - 1.0/n),  swr_(swr), snr_(snr)
        , betaNw_(betaNw), betaGn_(betaGn), betaGw_(betaGw), regardSnr_(regardSnr)
        {}

        Scalar alpha() const { return alpha_; }
        void setAlpha(Scalar alpha) { alpha_ = alpha; }

        Scalar m() const { return m_; }
        void setM(Scalar m) { m_ = m; n_ = 1.0/(1.0 - m); }

        Scalar n() const{ return n_; }
        void setN(Scalar n){ n_ = n; m_ = 1.0 - 1.0/n; }

        Scalar swr() const { return swr_; }
        void setSwr(Scalar swr) { swr_ = swr; }

        Scalar snr() const { return snr_; }
        void setSnr(Scalar swr) { snr_ = swr; }

        Scalar betaNw() const { return betaNw_; }
        void setBetaNw(Scalar betaNw) { betaNw_ = betaNw; }

        Scalar betaGn() const { return betaGn_; }
        void setBetaGn(Scalar betaGn) { betaGn_ = betaGn; }

        Scalar betaGw() const { return betaGw_; }
        void setBetaGw(Scalar betaGw) { betaGw_ = betaGw; }

        bool regardSnrForKrn() const { return regardSnr_; }
        void setRegardSnrForKrn(bool v) {regardSnr_ = v; }

        bool operator== (const Params& p) const
        {
            return Dune::FloatCmp::eq(alpha_, p.alpha_, 1e-6)
                   && Dune::FloatCmp::eq(n_, p.n_, 1e-6)
                   && Dune::FloatCmp::eq(m_, p.m_, 1e-6)
                   && Dune::FloatCmp::eq(swr_, p.swr_, 1e-6)
                   && Dune::FloatCmp::eq(snr_, p.snr_, 1e-6)
                   && Dune::FloatCmp::eq(betaNw_, p.betaNw_, 1e-6)
                   && Dune::FloatCmp::eq(betaGn_, p.betaGn_, 1e-6)
                   && Dune::FloatCmp::eq(betaGw_, p.betaGw_, 1e-6)
                   && regardSnr_ == p.regardSnr_;
        }

    private:
        Scalar alpha_, n_, m_, swr_, snr_, betaNw_, betaGn_, betaGw_;
        bool regardSnr_;
    };

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    template<class Scalar = double>
    static Params<Scalar> makeParams(const std::string& paramGroup)
    {
        const auto n = getParamFromGroup<Scalar>(paramGroup, "ParkerVanGenuchtenN");
        const auto alpha = getParamFromGroup<Scalar>(paramGroup, "ParkerVanGenuchtenAlpha");
        const auto swr = getParamFromGroup<Scalar>(paramGroup, "Swr", 0.0);
        const auto snr = getParamFromGroup<Scalar>(paramGroup, "Snr", 0.0);
        const auto betaNw = getParamFromGroup<Scalar>(paramGroup, "ParkerVanGenuchtenBetaNw", 1.0);
        const auto betaGn = getParamFromGroup<Scalar>(paramGroup, "ParkerVanGenuchtenBetaGn", 1.0);
        const auto betaGw = getParamFromGroup<Scalar>(paramGroup, "ParkerVanGenuchtenBetaGw", 1.0);
        const auto regardSnr = getParamFromGroup<bool>(paramGroup, "ParkerVanGenuchtenRegardSnrForKrn", false);
        return Params<Scalar>(alpha, n, swr, snr,
                              betaNw, betaGn, betaGw, regardSnr );
    }

    /*!
     * \brief The capillary pressure-saturation curve for the gas and wetting phase
     * \param params Array of parameters
     * \param swe Effective wetting phase saturation
     */
    template<class Scalar>
    static Scalar pcgw(Scalar swe, const Params<Scalar>& params)
    {
        assert(0 <= swe && swe <= 1);
        return pc_(swe, params);
    }

    /*!
     * \brief The capillary pressure-saturation curve for the non-wettigng and wetting phase
     * \param params Array of parameters
     * \param swe Effective wetting phase saturation
     */
    template<class Scalar>
    static Scalar pcnw(Scalar swe, const Params<Scalar>& params)
    {
        assert(0 <= swe && swe <= 1);
        return pc_(swe, params)/params.betaNw();
    }

    /*!
     * \brief The capillary pressure-saturation curve for the gas and nonwetting phase
     * \param params Array of parameters
     * \param ste Effective total liquid (wetting + nonwetting) saturation
     */
    template<class Scalar>
    static Scalar pcgn(const Scalar ste, const Params<Scalar>& params)
    {
        assert(0 <= ste && ste <= 1);
        return pc_(ste, params)/params.betaGn();
    }

    /*!
     * \brief This function ensures a continuous transition from 2 to 3 phases and vice versa
     * \param params Array of parameters
     * \param sne Non-wetting liquid saturation
     */
    template<class Scalar>
    static Scalar pcAlpha(Scalar sne, const Params<Scalar>& params)
    {
        /* regularization */
        if (sne <= 0.001)
            sne = 0.0;
        if (sne >= 1.0)
            sne = 1.0;

        if (sne > params.snr())
            return 1.0;
        else
        {
            if (params.snr() >= 0.001)
                return sne/params.snr();
            else
                return 0.0;
        };
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param swe Effective wetting phase saturation
     * \param params Array of parameters
     */
    template<class Scalar>
    static Scalar dpcgw_dswe(const Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        const Scalar powSeRegu = pow(swe, -1/params.m());
        return - 1.0/params.alpha() * pow(powSeRegu - 1, 1.0/params.n() - 1)/params.n()
            * powSeRegu/swe/params.m()/params.betaGw();
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param swe Effective wetting phase saturation
     * \param params Array of parameters
     */
    template<class Scalar>
    static Scalar dpcnw_dswe(const Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        const Scalar powSeRegu = pow(swe, -1/params.m());
        return - 1.0/params.alpha() * pow(powSeRegu - 1, 1.0/params.n() - 1)/params.n()
            * powSeRegu/swe/params.m()/params.betaNw();
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param ste Effective total liquid (wetting + nonwetting) saturation
     * \param params Array of parameters
     */
    template<class Scalar>
    static Scalar dpcgn_dste(const Scalar ste, const Params<Scalar>& params)
    {
        using std::pow;
        const Scalar powSeRegu = pow(ste, -1/params.m());
        return - 1.0/params.alpha() * pow(powSeRegu - 1, 1.0/params.n() - 1)/params.n()
            * powSeRegu/ste/params.m()/params.betaGn();
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
     * \param swe Effective wetting phase saturation
     * \param params Array of parameters.
     */
    template<class Scalar>
    static Scalar krw(const Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        using std::sqrt;
        const Scalar r = 1.0 - pow(1 - pow(swe, 1/params.m()), params.m());
        return sqrt(swe)*r*r;
    }

    /*!
     * \brief The relative permeability for the nonwetting phase
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
     * \param sn Absolute nonwetting liquid saturation
     * \param ste Effective total liquid (wetting + nonwetting) saturation
     */
    template<class Scalar>
    static Scalar krn(const Scalar swe, const Scalar sn, const Scalar ste, const Params<Scalar>& params)
    {
        Scalar krn;
        using std::pow;
        krn = pow(1 - pow(swe, 1/params.m()), params.m());
        krn -= pow(1 - pow(ste, 1/params.m()), params.m());
        krn *= krn;

        using std::clamp;
        using std::sqrt;
        if (params.regardSnrForKrn())
        {
            // regard Snr in the permeability of the n-phase, see Helmig1997
            const Scalar resIncluded = clamp(sn - params.snr()/ (1-params.swr()), 0.0, 1.0);
            krn *= sqrt(resIncluded);
        }
        else
            krn *= sqrt(sn / (1 - params.swr()));   // Hint: (ste - swe) = sn / (1-Swr)

        return krn;
    }

    /*!
     * \brief The relative permeability for the nonwetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * The permeability of gas in a 3p system equals the standard 2p description.
     * (see p61. in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83.) \cite delshad1989 <BR>
     *
     * \param params Array of parameters.
     * \param ste Effective total liquid (wetting + nonwetting) saturation
     */
    template<class Scalar>
    static Scalar krg(const Scalar ste, const Params<Scalar>& params)
    {
        assert(0 <= ste && ste <= 1);
        using std::cbrt;
        using std::pow;
        return cbrt(1 - ste) * pow(1 - pow(ste, 1/params.m()), 2*params.m());
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        gas phase in regard to the total liquid saturation of
     *        the medium as implied by the van Genuchten
     *        parameterization.
     *
     * \param ste The mobile total liquid saturation.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar dkrg_dste(const Scalar ste, const Params<Scalar>& params)
    {
        assert(0 < ste && ste <= 1);

        using std::pow;
        const Scalar x = pow(ste, 1.0/params.m());
        return
            -pow(1.0 - x, 2*params.m())
            *pow(1.0 - ste, -2.0/3)
            *(1.0/3 + 2*x/ste);
    }

    /*!
     * \brief The relative permeability for a phase.
     * \param params Array of parameters.
     * \param phaseIdx Indicator, The saturation of all phases.
     * \param swe Effective wetting phase saturation
     * \param sne Effective nonwetting saturation
     */
    template<class Scalar>
    static Scalar kr(const int phaseIdx, const Scalar swe, const Scalar sne, const Params<Scalar>& params)
    {
        switch (phaseIdx)
        {
            case 0:
                return krw(params, swe, sne);
            case 1:
                return krn(params, swe, sne);
            case 2:
                return krg(params, swe, sne);
        }
        DUNE_THROW(Dune::InvalidStateException,
                   "Invalid phase index ");
    }

private:

    /*!
     * \brief The standard van Genuchten two-phase pc-S relation either with respect to
     *        the effective wetting phase saturation Swe or the effective total liquid saturation Ste.
     * \param se Effective wetting phase or total liquid saturation
     * \param params Array of parameters.
     */
    template<class Scalar>
    const static Scalar pc_(const Scalar se, const Params<Scalar>& params)
    {
        using std::pow;
        return pow(pow(se, -1/params.m()) - 1, 1/params.n())/params.alpha();
    }

};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A regularization for the ParkerVanGenuchten3PRegularization material law
 * \note Regularization can be turned of by setting the threshold parameters
 *       out of range (runtime) or by replacing
 *       this class by NoRegularization (compile time).
 */
template <class Scalar>
class ParkerVanGenuchten3PRegularization
{
    using BaseLawParams = typename ParkerVanGenuchten3P::Params<Scalar>;

public:
    //! Regularization parameters
    template<class S>
    struct Params
    {
        /*!
         * \brief Set the threshold saturation below which the capillary pressure is regularized.
         *
         * Most problems are very sensitive to this value (e.g. making it smaller might
         * result in very high capillary pressures)
         */
        void setPcLowSwe(Scalar pcLowSwe)
        { pcLowSwe_ = pcLowSwe; }

        /*!
         * \brief Threshold saturation below which the capillary pressure is regularized.
         */
        Scalar pcLowSwe() const
        { return pcLowSwe_; }

        /*!
         * \brief Set the threshold saturation above which the capillary pressure is regularized.
         */
        void setPcHighSwe(Scalar pcHighSwe)
        { pcHighSwe_ = pcHighSwe; }

        /*!
         * \brief Threshold saturation above which the capillary pressure is regularized.
         *
         * Most problems are very sensitive to this value (e.g. making it smaller might
         * result in negative capillary pressures).
         */
        Scalar pcHighSwe() const
        { return pcHighSwe_; }

        /*!
         * \brief Set the threshold saturation below which the relative
         *        permeability of the nonwetting phase gets regularized.
         */
        void setKrnLowSwe(Scalar krnLowSwe)
        { krnLowSwe_ = krnLowSwe; }

        /*!
         * \brief Threshold saturation below which the relative
         *        permeability of the nonwetting phase gets regularized.
         */
        Scalar krnLowSwe() const
        { return krnLowSwe_; }

        /*!
         * \brief Set the threshold saturation below which the relative
         *        permeability of the nonwetting phase gets regularized.
         */
        void setKrgLowSte(Scalar krgLowSte)
        { krgLowSte_ = krgLowSte; }

        /*!
         * \brief Threshold saturation below which the relative
         *        permeability of the nonwetting phase gets regularized.
         */
        Scalar krgLowSte() const
        { return krgLowSte_; }

        /*!
         * \brief Set the threshold saturation above which the relative
         *        permeability of the wetting phase gets regularized.
         */
        void setKrwHighSwe(Scalar krwHighSwe)
        { krwHighSwe_ = krwHighSwe; }

        /*!
         * \brief Threshold saturation above which the relative
         *        permeability of the wetting phase gets regularized.
         */
        Scalar krwHighSwe() const
        { return krwHighSwe_; }


        /*!
         * \brief Choose whether to use a constant value for regularization of the
         *        pc-S curves or not
         * \param input True or false
         */
        void setConstRegularization(const bool input)
        { constRegularization_ = input; }

        /*!
         * \brief Returns whether to use a constant value for regularization of the
         *        pc-S curves or not
         */
        bool constRegularization() const
        { return constRegularization_; }

private:
        S pcLowSwe_ = 0.01;
        S pcHighSwe_ = 0.99;
        S krnLowSwe_ = 0.1;
        S krwHighSwe_ = 0.9;
        S krgLowSte_ = 1e-3;
        bool constRegularization_ = false;
    };

    //! Initialize the spline
    template<class MaterialLaw>
    void init(const MaterialLaw* m, const std::string& paramGroup)
    {
        pcLowSwe_ = getParamFromGroup<Scalar>(paramGroup, "ParkerVanGenuchtenPcLowSweThreshold", 0.01);
        pcHighSwe_ = getParamFromGroup<Scalar>(paramGroup, "ParkerVanGenuchtenPcHighSweThreshold", 0.99);
        krwHighSwe_ = getParamFromGroup<Scalar>(paramGroup, "ParkerVanGenuchtenKrwHighSweThreshold", 0.9);
        krnLowSwe_ = getParamFromGroup<Scalar>(paramGroup, "ParkerVanGenuchtenKrnLowSweThreshold", 0.1);
        krgLowSte_ = getParamFromGroup<Scalar>(paramGroup, "ParkerVanGenuchtenKrgLowSteThreshold", 1e-3);
        constRegularization_ = getParamFromGroup<bool>(paramGroup, "VanGenuchtenConstantRegularization", false);

        initPcParameters_(m, pcLowSwe_, pcHighSwe_);
        initKrParameters_(m, krnLowSwe_, krwHighSwe_);
    }

    template<class MaterialLaw, class BaseParams, class EffToAbsParams>
    void init(const MaterialLaw* m, const BaseParams& bp, const EffToAbsParams& etap, const Params<Scalar>& p)
    {
        pcLowSwe_ = p.pcLowSwe();
        pcHighSwe_ = p.pcHighSwe();
        krwHighSwe_ = p.krwHighSwe();
        krnLowSwe_ = p.krnLowSwe();
        krgLowSte_ = p.krgLowSte();
        constRegularization_ = p.constRegularization();

        initPcParameters_(m, pcLowSwe_, pcHighSwe_);
        initKrParameters_(m, krnLowSwe_, krwHighSwe_);
    }

    /*!
     * \brief Equality comparison with another instance
     */
    bool operator== (const ParkerVanGenuchten3PRegularization& o) const
    {
        return Dune::FloatCmp::eq(pcLowSwe_, o.pcLowSwe_, 1e-6)
               && Dune::FloatCmp::eq(pcHighSwe_, o.pcHighSwe_, 1e-6)
               && Dune::FloatCmp::eq(krwHighSwe_, o.krwHighSwe_, 1e-6)
               && Dune::FloatCmp::eq(krnLowSwe_, o.krnLowSwe_, 1e-6)
               && constRegularization_ == o.constRegularization_;
    }

    /*!
     * \brief The regularized capillary pressure-saturation curve for the gas and wetting phase
     * regularized part:
     *  - low saturation:  extend the \f$\mathrm{p_{cgw}(S_{we})}\f$ curve with the slope at the regularization point (i.e. no kink).
     *  - high saturation: connect the high regularization point with
     *                     with a spline and continue linearly for \f$\mathrm{S_{we} > 1}\f$
     * \param swe Effective wetting phase saturation
     */
    OptionalScalar<Scalar> pcgw(Scalar swe) const
    {
        // if specified, a constant value is used for regularization
        using std::clamp;
        if (constRegularization_)
            swe = clamp(swe, 0.0, 1.0);

        // make sure that the capillary pressure observes a derivative
        // != 0 for 'illegal' saturations. This is favourable for the
        // newton solver (if the derivative is calculated numerically)
        // in order to get the saturation moving to the right
        // direction if it temporarily is in an 'illegal' range.
        if (swe < pcLowSwe_)
            return pcgwLowSwePcgwValue_ + pcgwDerivativeLowSw_*(swe - pcLowSwe_);

        else if (swe > 1.0)
            return pcgwDerivativeHighSweEnd_*(swe - 1.0);

        else if (swe > pcHighSwe_)
            return pcgwSpline_.eval(swe);

        else
            return {}; // no regularization
    }

    /*!
     * \brief The regularized capillary pressure-saturation curve for the nonwetting and wetting phase
     * regularized part:
     *  - low saturation:  extend the \f$\mathrm{p_{cnw}(S_{we})}\f$ curve with the slope at the regularization point (i.e. no kink).
     *  - high saturation: connect the high regularization point with
     *                     with a spline and continue linearly for \f$\mathrm{S_{we} > 1}\f$
     * \param swe Effective wetting phase saturation
     */
    OptionalScalar<Scalar> pcnw(Scalar swe) const
    {
        // if specified, a constant value is used for regularization
        using std::clamp;
        if (constRegularization_)
            swe = clamp(swe, 0.0, 1.0);

        // make sure that the capillary pressure observes a derivative
        // != 0 for 'illegal' saturations. This is favourable for the
        // newton solver (if the derivative is calculated numerically)
        // in order to get the saturation moving to the right
        // direction if it temporarily is in an 'illegal' range.
        if (swe < pcLowSwe_)
            return pcnwLowSwePcnwValue_ + pcnwDerivativeLowSw_*(swe - pcLowSwe_);

         else if (swe > 1.0)
            return pcnwDerivativeHighSweEnd_*(swe - 1.0);

        else if (swe > pcHighSwe_)
            return pcnwSpline_.eval(swe);

        else
            return {}; // no regularization
    }

    /*!
     * \brief The regularized capillary pressure-saturation curve for the gas and nonwetting phase
     * regularized part:
     *  - low saturation:  extend the \f$\mathrm{p_{cgn}(S_{te})}\f$ curve with the slope at the regularization point (i.e. no kink).
     *  - high saturation: connect the high regularization point with
     *                     with a spline and continue linearly for \f$\mathrm{S_{te} > 1}\f$
     * \param ste Effective total liquid (sw + sn) saturation
     */
    OptionalScalar<Scalar> pcgn(Scalar ste) const
    {
        // if specified, a constant value is used for regularization
        using std::clamp;
        if (constRegularization_)
            ste = clamp(ste, 0.0, 1.0);


        // make sure that the capillary pressure observes a derivative
        // != 0 for 'illegal' saturations. This is favourable for the
        // newton solver (if the derivative is calculated numerically)
        // in order to get the saturation moving to the right
        // direction if it temporarily is in an 'illegal' range.
        const Scalar pcLowSte = pcLowSwe_;
        const Scalar pcHighSte = pcHighSwe_;
        if (ste < pcLowSte)
            return pcgnLowStePcgnValue_ + pcgnDerivativeLowSt_*(ste - pcLowSte);

        else if (ste > 1.0)
            return pcgnDerivativeHighSteEnd_*(ste - 1.0);

        else if (ste > pcHighSte)
            return pcgnSpline_.eval(ste);

        else
            return {}; // no regularization
    }

    /*!
     * \brief This function ensures a continuous transition from 2 to 3 phases and vice versa
     * \param sne Effective nonwetting liquid saturation
     */
    OptionalScalar<Scalar> pcAlpha(Scalar sne) const
    {
        // no regularization
        return {};
    }

    /*!
     * \brief The regularized relative permeability for the wetting phase
     * \param swe Effective wetting phase saturation
     */
    OptionalScalar<Scalar> krw(const Scalar swe) const
    {
        if (swe < 0.0)
            return 0.0;
        else if (swe > 1.0 - std::numeric_limits<Scalar>::epsilon())
            return 1.0;
        else
            return {}; // no regularization
    }


    /*!
     * \brief The regularized relative permeability for the nonwetting phase
     * \param swe Effective wetting phase saturation
     * \param sn Nonwetting saturation
     * \param ste Effective total (wetting + nonwetting) saturation
     */
    OptionalScalar<Scalar> krn(Scalar swe, const Scalar sn, Scalar ste) const
    {
        using std::clamp;
        swe = clamp(swe, 0.0, 1.0);
        ste = clamp(ste, 0.0, 1.0);

        if (ste - swe <= 0.0)
            return 0.0;
        else
            return ParkerVanGenuchten3P::krn(swe, sn, ste, *baseLawParamsPtr_);
    }

    /*!
     * \brief The regularized relative permeability for the gas phase
     * \param ste Effective total (wetting + nonwetting) saturation
     */
    OptionalScalar<Scalar> krg(const Scalar ste) const
    {
        //return 0 if there is no gas
        if (ste > 1.0 - std::numeric_limits<Scalar>::epsilon())
            return 0.0;

        // use linear regularization for very high gas saturations
        // to avoid a kink in the curve and to maintain a slope for
        // the Newton solver
        if (ste <= krgLowSte_)
            return krgLowStkrgValue_ + krgDerivativeLowSt_*(ste - krgLowSte_);
        else
        {
            // For very low gas saturations:
            // We use a scaling factor that decreases the gas phase permeability quite fast a very low gas phase
            // saturations, thus making that phase virtually immobile.
            // This prevents numerical issues related to the degeneration of the gas phase mass balance for the 3p3c model
            // at very low gas phase saturations.

            // get the absolute gas phase saturation
            const Scalar st = ste*(1 - swr_) + swr_;
            const Scalar sg = 1.0 - st;

            // do not regularize
            if (sg > 0.1)
                return {};

            // return original curve scaled by factor
            using std::max;
            const Scalar scalFact = max(0.0, (sg - sgr_)/(0.1 - sgr_));

            return ParkerVanGenuchten3P::krg(ste, *baseLawParamsPtr_) * scalFact;
        }
    }

    /*!
     * \brief The relative permeability for a phase.
     * \param phaseIdx Indicator, The saturation of all phases.
     * \param swe Effective wetting phase saturation
     * \param sne Effective nonwetting saturation
     */
    OptionalScalar<Scalar> kr(const int phaseIdx, const Scalar swe, const Scalar sne) const
    {
        switch (phaseIdx)
        {
        case 0:
            return krw(swe, sne);
        case 1:
            return krn(swe, sne);
        case 2:
            return krg(swe, sne);
        }
        DUNE_THROW(Dune::InvalidStateException,
                   "Invalid phase index ");
    }

private:
    template<class MaterialLaw>
    void initPcParameters_(const MaterialLaw* m, const Scalar lowSwe, const Scalar highSwe)
    {
        const auto lowSw = MaterialLaw::EffToAbs::sweToSw(lowSwe, m->effToAbsParams());
        const auto highSw = MaterialLaw::EffToAbs::sweToSw(highSwe, m->effToAbsParams());
        const auto dsw_dswe = MaterialLaw::EffToAbs::dsw_dswe(m->effToAbsParams());
        const auto dst_dste = MaterialLaw::EffToAbs::dst_dste(m->effToAbsParams());

        baseLawParamsPtr_ = &m->basicParams();

        // pcgw
        pcgwLowSwePcgwValue_ = m->template pcgw<false>(lowSw, 0.0);
        pcgwDerivativeLowSw_ = m->template dpcgw_dsw<false>(lowSw, 0.0)*dsw_dswe;
        pcgwHighSwePcgwValue_ = m->template pcgw<false>(highSw, 0.0);
        pcgwDerivativeHighSweThreshold_ = m->template dpcgw_dsw<false>(highSw, 0.0)*dsw_dswe;
        pcgwDerivativeHighSweEnd_ = 2.0*(0.0 - m->template pcgw<false>(highSw, 0.0))/(1.0 - highSwe);
        pcgwSpline_ = Spline<Scalar>(highSwe, 1.0, // x0, x1
                                     pcgwHighSwePcgwValue_, 0, // y0, y1
                                     pcgwDerivativeHighSweThreshold_, pcgwDerivativeHighSweEnd_); // m0, m1

        // pcnw
        pcnwLowSwePcnwValue_ = m->template pcnw<false>(lowSw, 0.0);
        pcnwDerivativeLowSw_ = m->template dpcnw_dsw<false>(lowSw, 0.0)*dsw_dswe;
        pcnwHighSwePcnwValue_ = m->template pcnw<false>(highSw, 0.0);
        pcnwDerivativeHighSweThreshold_ = m->template dpcnw_dsw<false>(highSw, 0.0);
        pcnwDerivativeHighSweEnd_ = 2.0*(0.0 - m->template pcnw<false>(highSw, 0.0))/(1.0 - highSwe);
        pcnwSpline_ = Spline<Scalar>(highSwe, 1.0, // x0, x1
                                     pcnwHighSwePcnwValue_, 0, // y0, y1
                                     pcnwDerivativeHighSweThreshold_, pcnwDerivativeHighSweEnd_); // m0, m1

        // pcgn
        pcgnLowStePcgnValue_ = m->template pcgn<false>(lowSw, 0.0);
        pcgnDerivativeLowSt_ = m->template dpcgn_dst<false>(lowSw, 0.0)*dst_dste;
        pcgnHighSwePcgnValue_ = m->template pcgn<false>(highSw, 0.0);
        pcgnDerivativeHighSteThreshold_ = m->template dpcgn_dst<false>(highSw, 0.0);
        pcgnDerivativeHighSteEnd_ = 2.0*(0.0 - m->template pcgn<false>(highSw, 0.0))/(1.0 - highSwe);
        pcgnSpline_ = Spline<Scalar>(highSwe, 1.0, // x0, x1
                                     pcgnHighSwePcgnValue_, 0, // y0, y1
                                     pcgnDerivativeHighSteThreshold_, pcgnDerivativeHighSteEnd_); // m0, m1

    }

    template<class MaterialLaw>
    void initKrParameters_(const MaterialLaw* m, const Scalar lowSwe, const Scalar highSwe)
    {
        krgLowStkrgValue_ = ParkerVanGenuchten3P::krg(krgLowSte_, *baseLawParamsPtr_);
        krgDerivativeLowSt_ = ParkerVanGenuchten3P::dkrg_dste(krgLowSte_, *baseLawParamsPtr_);

        swr_ = m->effToAbsParams().swr();
        sgr_ = m->effToAbsParams().sgr();
    }

    Scalar krgLowStkrgValue_;
    Scalar krgDerivativeLowSt_;

    Scalar pcLowSwe_, pcHighSwe_;
    Scalar krwHighSwe_, krnLowSwe_, krgLowSte_;

    // pcgw
    Scalar pcgwLowSwePcgwValue_;
    Scalar pcgwHighSwePcgwValue_;
    Scalar pcgwDerivativeLowSw_;
    Scalar pcgwDerivativeHighSweThreshold_;
    Scalar pcgwDerivativeHighSweEnd_;

    // pcgn
    Scalar pcgnLowStePcgnValue_;
    Scalar pcgnHighSwePcgnValue_;
    Scalar pcgnDerivativeLowSt_;
    Scalar pcgnDerivativeHighSteThreshold_;
    Scalar pcgnDerivativeHighSteEnd_;

    // pcnw
    Scalar pcnwLowSwePcnwValue_;
    Scalar pcnwHighSwePcnwValue_;
    Scalar pcnwDerivativeLowSw_;
    Scalar pcnwDerivativeHighSweThreshold_;
    Scalar pcnwDerivativeHighSweEnd_;

    Spline<Scalar> pcgwSpline_;
    Spline<Scalar> pcnwSpline_;
    Spline<Scalar> pcgnSpline_;
    Spline<Scalar> krwSpline_;
    Spline<Scalar> krnSpline_;

    Scalar swr_, sgr_;

    bool constRegularization_;

    const BaseLawParams* baseLawParamsPtr_;
};

template<class ScalarType,
         class BaseLaw,
         class Regularization = NoRegularization,
         class EffToAbsPolicy = ParkerVanGenuchten3PEffToAbsPolicy>
class ParkerVanGenuchtenMaterialLaw : public Adapter<ParkerVanGenuchtenMaterialLaw<ScalarType, BaseLaw, Regularization, EffToAbsPolicy>, ThreePhasePcKrSw>
{
public:

    using Scalar = ScalarType;

    using BasicParams = typename BaseLaw::template Params<Scalar>;
    using EffToAbsParams = typename EffToAbsPolicy::template Params<Scalar>;
    using RegularizationParams = typename Regularization::template Params<Scalar>;

    using EffToAbs = EffToAbsPolicy;

    /*!
     * \brief Return whether this law is regularized
     */
    static constexpr bool isRegularized()
    { return !std::is_same<Regularization, NoRegularization>::value; }

    /*!
     * \brief Deleted default constructor (so we are never in an undefined state)
     * \note store owning pointers to laws instead if you need default-constructible objects
     */
    ParkerVanGenuchtenMaterialLaw() = delete;

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    explicit ParkerVanGenuchtenMaterialLaw(const std::string& paramGroup)
    : basicParams_(makeBasicParams(paramGroup))
    , effToAbsParams_(makeEffToAbsParams(paramGroup))
    {
        if constexpr (isRegularized())
            regularization_.init(this, paramGroup);
    }

    /*!
     * \brief Construct from parameter structs
     * \note More efficient constructor but you need to ensure all parameters are initialized
     */
    ParkerVanGenuchtenMaterialLaw(const BasicParams& baseParams,
                                  const EffToAbsParams& effToAbsParams = {},
                                  const RegularizationParams& regParams = {})
    : basicParams_(baseParams)
    , effToAbsParams_(effToAbsParams)
    {
        if constexpr (isRegularized())
            regularization_.init(this, baseParams, effToAbsParams, regParams);
    }

    /*!
     * \brief The capillary pressure-saturation curve for the gas and wetting phase
     */
    template<bool enableRegularization = isRegularized()>
    Scalar pcgw(const Scalar sw, const Scalar /*dummySn*/) const
    {
        const auto swe = EffToAbs::swToSwe(sw, effToAbsParams_);
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.pcgw(swe);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::pcgw(swe, basicParams_);
    }

    /*!
     * \brief The capillary pressure-saturation curve for the nonwetting and wetting phase
     */
    template<bool enableRegularization = isRegularized()>
    Scalar pcnw(const Scalar sw, const Scalar /*dummySn*/) const
    {
        const auto swe = EffToAbs::swToSwe(sw, effToAbsParams_);
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.pcnw(swe);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::pcnw(swe, basicParams_);
    }

    /*!
     * \brief The capillary pressure-saturation curve for the gas and nonwetting phase
     * \param sw Wetting saturation
     * \param sn Nonwetting saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar pcgn(const Scalar sw, const Scalar sn) const
    {
        const auto swe = EffToAbs::swToSwe(sw + sn, effToAbsParams_);
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.pcgn(swe);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::pcgn(swe, basicParams_);
    }

    /*!
     * \brief This function ensures a continuous transition from 2 to 3 phases and vice versa
     */
    template<bool enableRegularization = isRegularized()>
    Scalar pcAlpha(const Scalar /*dummySw*/, const Scalar sn) const
    {
        const auto sne = EffToAbs::snToSne(sn, effToAbsParams_);
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.pcAlpha(sne);
            if (regularized)
                return regularized.value();
        }
        return BaseLaw::pcAlpha(sne, basicParams_);
    }

    /*!
     * \brief The partial derivative of the capillary pressure w.r.t. the saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar dpcgw_dsw(const Scalar sw, const Scalar /*dummyS*/) const
    {
        const auto swe = EffToAbs::swToSwe(sw, effToAbsParams_);
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.dpcgw_dswe(swe);
            if (regularized)
                return regularized.value()*EffToAbs::dswe_dsw(effToAbsParams_);
        }

        return BaseLaw::dpcgw_dswe(swe, basicParams_)*EffToAbs::dswe_dsw(effToAbsParams_);
    }

    /*!
     * \brief The partial derivative of the capillary pressure w.r.t. the saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar dpcnw_dsw(const Scalar sw, const Scalar /*dummySw*/) const
    {
        const auto swe = EffToAbs::swToSwe(sw, effToAbsParams_);
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.dpcnw_dswe(swe);
            if (regularized)
                return regularized.value()*EffToAbs::dswe_dsw(effToAbsParams_);
        }

        return BaseLaw::dpcnw_dswe(swe, basicParams_)*EffToAbs::dswe_dsw(effToAbsParams_);
    }

    /*!
     * \brief The partial derivative of the capillary pressure w.r.t. the saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar dpcgn_dst(const Scalar st, const Scalar /*dummySw*/) const
    {
        const auto ste = EffToAbs::stToSte(st, effToAbsParams_);
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.dpcgn_dste(ste);
            if (regularized)
                return regularized.value()*EffToAbs::dswte_dst(effToAbsParams_);
        }

        return BaseLaw::dpcgn_dste(ste, basicParams_)*EffToAbs::dste_dst(effToAbsParams_);
    }

    /*!
     * \brief The relative permeability for the wetting phase
     * \param sw Wetting saturation
     * \param sn Nonwetting saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar krw(const Scalar sw, const Scalar sn) const
    {
        const auto swe = EffToAbs::swToSwe(sw, effToAbsParams_);
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.krw(swe);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::krw(swe, basicParams_);
    }

    /*!
     * \brief The relative permeability for the nonwetting phase
     * \param sw Wetting saturation
     * \param sn Nonwetting saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar krn(const Scalar sw, const Scalar sn) const
    {
        const auto swe = EffToAbs::swToSwe(sw, effToAbsParams_);
        const auto ste = EffToAbs::stToSte(sw + sn, effToAbsParams_);
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.krn(swe, sn, ste);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::krn(swe, sn, ste, basicParams_);
    }

    /*!
     * \brief The relative permeability for the nonwetting phase
     * \param sw Wetting saturation
     * \param sn Nonwetting saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar krg(const Scalar sw, const Scalar sn) const
    {
        const auto ste = EffToAbs::stToSte(sw + sn, effToAbsParams_);
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.krg(ste);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::krg(ste, basicParams_);
    }

    /*!
     * \brief The relative permeability for the nonwetting phase
     * \param phaseIdx Indicator, The saturation of all phases.
     * \param sw Wetting saturation
     * \param sn Nonwetting saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar kr(const int phaseIdx, const Scalar sw, const Scalar sn) const
    {
        switch (phaseIdx)
        {
            case 0:
                return krw(sw, sn);
            case 1:
                return krn(sw, sn);
            case 2:
                return krg(sw, sn);
        }
        DUNE_THROW(Dune::InvalidStateException,
                   "Invalid phase index ");
    }

    /*!
     * \brief The derivative of the relative permeability for the nonwetting phase w.r.t. saturation
     * \param st Total (wetting + nonwetting) saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar dkrg_dst(const Scalar st) const
    {
        const auto ste = EffToAbs::stToSte(st, effToAbsParams_);
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.dkrg_dste(ste);
            if (regularized)
                return regularized.value()*EffToAbs::dste_dst(effToAbsParams_);
        }

        return BaseLaw::dkrg_dste(ste, basicParams_)*EffToAbs::dste_dst(effToAbsParams_);
    }

    /*!
     * \brief Equality comparison with another instance
     */
    bool operator== (const ParkerVanGenuchtenMaterialLaw& o) const
    {
        return basicParams_ == o.basicParams_
               && effToAbsParams_ == o.effToAbsParams_
               && regularization_ == o.regularization_;
    }

    /*!
     * \brief Create the base law's parameters using
     *        input file parameters
     */
    static BasicParams makeBasicParams(const std::string& paramGroup)
    {
        return BaseLaw::template makeParams<Scalar>(paramGroup);
    }

    /*!
     * \brief Return the base law's parameters
     */
    const BasicParams& basicParams() const
    { return basicParams_; }

    /*!
     * \brief Create the parameters of the EffToAbs policy using
     *        input file parameters
     */
    static EffToAbsParams makeEffToAbsParams(const std::string& paramGroup)
    {
        return EffToAbs::template makeParams<Scalar>(paramGroup);
    }

    /*!
     * \brief Return the parameters of the EffToAbs policy
     */
    const EffToAbsParams& effToAbsParams() const
    { return effToAbsParams_; }

private:
    BasicParams basicParams_;
    EffToAbsParams effToAbsParams_;
    Regularization regularization_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A configuration for using the ParkerVanGenuchten material law without regularization
 */
template<class Scalar>
using ParkerVanGenuchten3PNoReg = ParkerVanGenuchtenMaterialLaw<Scalar, ParkerVanGenuchten3P, NoRegularization, ParkerVanGenuchten3PEffToAbsPolicy>;

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A default configuration for using the ParkerVanGenuchten material law
 */
template<class Scalar>
using ParkerVanGenuchten3PDefault = ParkerVanGenuchtenMaterialLaw<Scalar, ParkerVanGenuchten3P, ParkerVanGenuchten3PRegularization<Scalar>, ParkerVanGenuchten3PEffToAbsPolicy>;

} // end namespace Dumux::FluidMatrix

#endif
