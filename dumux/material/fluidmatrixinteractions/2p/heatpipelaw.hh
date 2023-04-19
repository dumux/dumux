// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the capillary pressure <-> saturation relation
 *        for the heatpipe problem.
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWOP_HEATPIPELAW_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_TWOP_HEATPIPELAW_HH

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
     * \param sw Saturation of the wetting phase \f$\mathrm{S_w}\f$
     */
    Scalar pc(const Scalar sw) const
    {
        const Scalar swe = EffToAbs::swToSwe(sw, effToAbsParams_);
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
    Scalar endPointPc() const
    { return 0.0; }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation.
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{S_w}\f$
     */
    Scalar dpc_dsw(const Scalar sw) const
    {
        const Scalar swe = EffToAbs::swToSwe(sw, effToAbsParams_);
        const Scalar sne = 1 - swe;
        const Scalar p0Gamma = params_.p0()*params_.gamma();
        if (sne > 1.0)
            sne = 1.0;
        else if (sne <= 0.0)
            return -p0Gamma*1.417*EffToAbs::dswe_dsw(effToAbsParams_);
        else
            return - p0Gamma*((3*1.263*sne - 2*2.120)*sne + 1.417)*EffToAbs::dswe_dsw(effToAbsParams_);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium.
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{S_w}\f$
     */
    Scalar krw(const Scalar sw) const
    {
        const Scalar swe = EffToAbs::swToSwe(sw, effToAbsParams_);
        return kr_(swe);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium.
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{S_w}\f$
     */
    Scalar krn(const Scalar sw) const
    {
        const Scalar swe = EffToAbs::swToSwe(sw, effToAbsParams_);
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
