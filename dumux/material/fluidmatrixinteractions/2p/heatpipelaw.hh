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
 * \brief Implementation of the capillary pressure <-> saturation relation
 *        for the heatpipe problem.
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWOP_HEATPIPELAW_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_TWOP_HEATPIPELAW_HH

// remove from here after release 3.3 /////////////
#include "heatpipelawparams.hh"

#include <dumux/common/spline.hh>

#include <algorithm>

#include <math.h>
#include <assert.h>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the capillary pressure <-> saturation
 *        relation for the heatpipe problem.
 *
 * This class bundles the "raw" curves as static members and doesn't concern itself
 * converting absolute to effective saturations and vince versa.
 */
template <class ScalarT, class ParamsT = HeatPipeLawParams<ScalarT> >
class [[deprecated("Use new material laws and FluidMatrix::HeatPipeLaw instead!")]] HeatPipeLaw
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     * \param params Array of parameters asd
     * \param Sw Effective saturation of of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     */
    static Scalar pc(const Params &params, Scalar Sw)
    {
        Scalar Sn = 1 - Sw;
        Scalar p0Gamma = params.p0()*params.gamma();

        // regularization
        if (Sn >= 1.0) {
            Scalar y = p0Gamma*(  (1.263*1.0 -   2.120)*1.0 + 1.417)*1.0;
            Scalar m = p0Gamma*((3*1.263*1.0 - 2*2.120)*1.0 + 1.417);
            return (Sn - 1)*m + y;
        }
        else if (Sn <= 0.0) {
            Scalar y = 0.0;
            Scalar m = p0Gamma*1.417;
            return Sn*m + y;
        }

        return p0Gamma*((1.263*Sn - 2.120)*Sn + 1.417) * Sn;
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * \return The effective saturaion of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params Array of parameters
     * \param pC Capillary pressure \f$\mathrm{[p_C]}\f$ in \f$\mathrm{[Pa]}\f$.
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        DUNE_THROW(Dune::NotImplemented, "HeatPipeLaw::Sw");
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param params Array of parameters
     * \param Sw Effective saturation of of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     */
    static Scalar dpC_dSw(const Params &params, Scalar Sw)
    {
        Scalar Sn = 1 - Sw;
        Scalar p0Gamma = params.p0()*params.gamma();
        if (Sn > 1.0)
            Sn = 1.0;
        else if (Sn <= 0.0) {
            Scalar m = -p0Gamma*1.417;
            return m;
        }

        Scalar m = - p0Gamma*((3*1.263*Sn - 2*2.120)*Sn + 1.417);
        return m;
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     * \param params Array of parameters
     * \param pC Capillary pressure \f$\mathrm{[p_C]}\f$ in \f$\mathrm{[Pa]}\f$.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        DUNE_THROW(Dune::NotImplemented, "HeatPipeLaw::dSw_dpC");
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param params Array of parameters
     * \param Sw The mobile saturation of the wetting phase.
     */
    static Scalar krw(const Params &params, Scalar Sw)
    {
        return kr_(Sw);
    }

    /*!
     * \brief The relative permeability for the nonwetting phase.
     *
     * \param params Array of parameters
     * \param Sw The mobile saturation of the wetting phase.
     */
    static Scalar krn(const Params &params, Scalar Sw)
    {
        Scalar Sn = 1 - Sw;
        return kr_(Sn);
    }

private:
    static Scalar kr_(Scalar S)
    {
        const Scalar eps = 0.95;
        if (S >= 1)
            return 1;
        else if (S <= 0)
            return 0;
        else if (S > eps) {
            // regularize
            using Spline = Dumux::Spline<Scalar>;
            Spline sp(eps, 1.0, // x1, x2
                      eps*eps*eps, 1, // y1, y2
                      3*eps*eps, 0); // m1, m2
            return sp.eval(S);
        }

        return S*S*S;
    }
};

} // end namespace Dumux
// remove until here after release 3.3 /////////////

#include <cmath>
#include <algorithm>

#include <dumux/common/parameters.hh>
#include <dumux/common/spline.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabsdefaultpolicy.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the capillary pressure <-> saturation
 *        relation for the heatpipe problem.
 */
template<class ScalarType, class EffToAbsPolicy = TwoPEffToAbsDefaultPolicy>
class HeatPipeLaw : Adapter<HeatPipeLaw<ScalarType, EffToAbsPolicy>, PcKrSw>
{

public:
    using Scalar = ScalarType;
    using EffToAbsParams = typename EffToAbsPolicy::template Params<Scalar>;
    using EffToAbs = EffToAbsPolicy;

    /*!
     * \brief Return the number of fluid phases
     */
    static constexpr int numFluidPhases()
    { return 2; }

    /*!
     * \brief Return whether this law is regularized
     */
    static constexpr bool isRegularized()
    { return false; }

    /*!
     * \brief The parameter type
     * \tparam Scalar The scalar type
     */
    struct Params
    {
        Params(Scalar gamma, Scalar p0)
        : gamma_(gamma), p0_(p0)
        {}

        Scalar gamma() const { return gamma_; }
        void setGamma(Scalar gamma) { gamma_ = gamma; }

        Scalar p0() const { return p0_; }
        void setP0(Scalar p0) { p0_ = p0; }

        bool operator== (const Params& p) const
        {
            return Dune::FloatCmp::eq(gamma(), p.gamma(), 1e-6)
                   && Dune::FloatCmp::eq(p0(), p.p0(), 1e-6);
        }

    private:
        Scalar gamma_, p0_;
    };

    /*!
     * \brief Deleted default constructor (so we are never in an undefined state)
     * \note store owning pointers to laws instead if you need default-constructible objects
     */
    HeatPipeLaw() = delete;

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    explicit HeatPipeLaw(const std::string& paramGroup)
    {
        params_ = makeParams(paramGroup);
        effToAbsParams_ = EffToAbs::template makeParams<Scalar>(paramGroup);
    }

    /*!
     * \brief Construct from parameter structs
     * \note More efficient constructor but you need to ensure all parameters are initialized
     */
    HeatPipeLaw(const Params& params,
                const EffToAbsParams& effToAbsParams = {})
    : params_(params)
    , effToAbsParams_(effToAbsParams)
    , spline_(eps_, 1.0, // x1, x2
              eps_*eps_*eps_, 1.0, // y1, y2
              3.0*eps_*eps_, 0.0) // m1, m2
    {}

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    static Params makeParams(const std::string& paramGroup)
    {
        const auto gamma = getParamFromGroup<Scalar>(paramGroup, "HeatpipeLawGamma");
        const auto p0 = getParamFromGroup<Scalar>(paramGroup, "HeatpipeLawP0");
        return {gamma, p0};
    }

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     * \param swe Effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     */
    template<class Scalar>
    Scalar pc(Scalar swe) const
    {
        const Scalar sne = 1 - swe;
        const Scalar p0Gamma = params_.p0()*params_.gamma();

        // regularization
        if (sne >= 1.0)
        {
            const Scalar y = p0Gamma*(  (1.263*1.0 -   2.120)*1.0 + 1.417)*1.0;
            const Scalar m = p0Gamma*((3*1.263*1.0 - 2*2.120)*1.0 + 1.417);
            return (sne - 1)*m + y;
        }
        else if (sne <= 0.0)
            return sne * p0Gamma*1.417;
        else
            return p0Gamma*((1.263*sne - 2.120)*sne + 1.417) * sne;
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     */
    template<class Scalar>
    Scalar endPointPc() const
    { return 0.0; }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation.
     *
     * \param swe Effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     */
    template<class Scalar>
    Scalar dpc_dswe(Scalar swe) const
    {
        const Scalar sne = 1 - swe;
        const Scalar p0Gamma = params_.p0()*params_.gamma();
        if (sne > 1.0)
            sne = 1.0;
        else if (sne <= 0.0)
            return -p0Gamma*1.417;
        else
            return - p0Gamma*((3*1.263*sne - 2*2.120)*sne + 1.417);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium.
     *
     * \param swe The mobile saturation of the wetting phase.
     */
    template<class Scalar>
    Scalar krw(Scalar swe) const
    {
        return kr_(swe);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium.
     *
     * \param swe The mobile saturation of the wetting phase.
     */
    template<class Scalar>
    Scalar krn(Scalar swe) const
    {
        const Scalar sne = 1 - swe; // TODO does this make sense?
        return kr_(sne);
    }

    /*!
     * \brief Equality comparison with another instance
     */
    bool operator== (const HeatPipeLaw<Scalar, EffToAbs>& o) const
    {
        return params_ == o.params_
               && effToAbsParams_ == o.effToAbsParams_;
    }

    const EffToAbsParams& effToAbsParams() const
    { return effToAbsParams_; }

private:

    Scalar kr_(Scalar s) const
    {
        if (s >= 1.0)
            return 1;
        else if (s <= 0.0)
            return 0;
        else if (s > eps_)
            return spline_.eval(s); // regularize
        else
            return s*s*s;
    }

    Params params_;
    EffToAbsParams effToAbsParams_;
    static constexpr Scalar eps_ = 0.95;
    Dumux::Spline<Scalar> spline_;
};
} // end namespace Dumux::FluidMatrix

#endif
