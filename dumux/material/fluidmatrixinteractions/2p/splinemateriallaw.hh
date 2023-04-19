// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief A spline approximation wrapper for 2p material laws
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWOP_SPLINE_MATERIAL_LAW_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_TWOP_SPLINE_MATERIAL_LAW_HH

#include <memory> // unique_ptr

#include <dumux/common/parameters.hh>
#include <dumux/common/monotonecubicspline.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/material/fluidmatrixinteractions/2p/materiallaw.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A spline approximation wrapper for 2p material laws
 * \tparam TwoPMaterialLaw the type of material law to be wrapped
 * \tparam approximatePcSwInverse if this is set true, the
 *          spline approximates sw(pc) and evaluating pc(sw) needs spline inversion.
 *          if this is false, the spline approximates pc(sw) and evaluating
 *          sw(pc) needs spline inversion. Spline inversion is rather expensive
 *          since it has to be done numerically.
 */
template<class TwoPMaterialLaw, bool approximatePcSwInverse = false>
class SplineTwoPMaterialLaw
: public TwoPMaterialLaw
, public Adapter<SplineTwoPMaterialLaw<TwoPMaterialLaw, approximatePcSwInverse>, PcKrSw>
{
public:
    using Scalar = typename TwoPMaterialLaw::Scalar;

    using BasicParams = typename TwoPMaterialLaw::BasicParams;
    using EffToAbsParams = typename TwoPMaterialLaw::BasicParams;
    using RegularizationParams = typename TwoPMaterialLaw::RegularizationParams;

    /*!
     * \brief Return the number of fluid phases
     */
    static constexpr int numFluidPhases()
    { return 2; }

    /*!
     * \brief We are always regularized in the sense that we replace
     *        the original curve by a cubic spline
     */
    static constexpr bool isRegularized()
    { return true; }

    /*!
     * \brief Deleted default constructor (so we are never in an undefined state)
     * \note store owning pointers to laws instead if you need default-constructible objects
     */
    SplineTwoPMaterialLaw() = delete;

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    explicit SplineTwoPMaterialLaw(const std::string& paramGroup)
    : TwoPMaterialLaw(paramGroup)
    {
        const std::array<Scalar, 2> defaultInterval{{ 0.01, 1.0 }};
        const auto sweInterval = getParamFromGroup<std::array<Scalar, 2>>(paramGroup, "SplineSweInterval", defaultInterval);
        swInterval_ = {{ TwoPMaterialLaw::EffToAbs::sweToSw(sweInterval[0], this->effToAbsParams()),
                         TwoPMaterialLaw::EffToAbs::sweToSw(sweInterval[1], this->effToAbsParams()) }};
        swIntervalPc_ = {{ TwoPMaterialLaw::pc(swInterval_[1]),
                           TwoPMaterialLaw::pc(swInterval_[0]) }};
        numSwSamples_ = getParamFromGroup<std::size_t>(paramGroup, "SplineNumSwSamples", 30);

        pcSpline_ = makeSweSpline_(
            [&](const Scalar s){ return TwoPMaterialLaw::pc(s); },
            approximatePcSwInverse
        );

        krwSpline_ = makeSweSpline_(
            [&](const Scalar s){ return TwoPMaterialLaw::krw(s); }
        );

        krnSpline_ = makeSweSpline_(
            [&](const Scalar s){ return TwoPMaterialLaw::krn(s); }
        );
    }

    /*!
     * \brief Construct from parameter structs
     * \note More efficient constructor but you need to ensure all parameters are initialized
     */
    SplineTwoPMaterialLaw(const std::array<Scalar, 2>& sweInterval,
                          std::size_t numSwSamples,
                          TwoPMaterialLaw&& twoP)
    : TwoPMaterialLaw(std::move(twoP))
    , numSwSamples_(numSwSamples)
    {
        swInterval_ = {{ TwoPMaterialLaw::EffToAbs::sweToSw(sweInterval[0], this->effToAbsParams()),
                         TwoPMaterialLaw::EffToAbs::sweToSw(sweInterval[1], this->effToAbsParams()) }};
        swIntervalPc_ = {{ TwoPMaterialLaw::pc(swInterval_[1]),
                           TwoPMaterialLaw::pc(swInterval_[0]) }};
    }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    Scalar pc(const Scalar sw) const
    {
        if (sw > swInterval_[0] && sw < swInterval_[1])
        {
            if constexpr (approximatePcSwInverse)
                return pcSpline_->evalInverse(sw);
            else
                return pcSpline_->eval(sw);
        }

        return TwoPMaterialLaw::pc(sw);
    }

    /*!
     * \brief The partial derivative of the capillary pressure w.r.t. the saturation
     */
    Scalar dpc_dsw(const Scalar sw) const
    {
        if (sw > swInterval_[0] && sw < swInterval_[1])
        {
            if constexpr (approximatePcSwInverse)
                return 1.0/pcSpline_->evalDerivative(pcSpline_->evalInverse(sw));
            else
                return pcSpline_->evalDerivative(sw);
        }

        return TwoPMaterialLaw::dpc_dsw(sw);
    }

    /*!
     * \brief The saturation-capillary pressure curve
     */
    Scalar sw(const Scalar pc) const
    {
        if (pc > swIntervalPc_[0] && pc < swIntervalPc_[1])
        {
            if constexpr (approximatePcSwInverse)
                return pcSpline_->eval(pc);
            else
                return pcSpline_->evalInverse(pc);
        }

        return TwoPMaterialLaw::sw(pc);
    }

    /*!
     * \brief The partial derivative of the saturation to the capillary pressure
     */
    Scalar dsw_dpc(const Scalar pc) const
    {
        if (pc > swIntervalPc_[0] && pc < swIntervalPc_[1])
        {
            if constexpr (approximatePcSwInverse)
                return pcSpline_->evalDerivative(pc);
            else
                return 1.0/pcSpline_->evalDerivative(pcSpline_->evalInverse(pc));
        }

        return TwoPMaterialLaw::dsw_dpc(pc);
    }

    /*!
     * \brief The relative permeability for the wetting phase
     */
    Scalar krw(const Scalar sw) const
    {
        if (sw > swInterval_[0] && sw < swInterval_[1])
            return krwSpline_->eval(sw);

        return TwoPMaterialLaw::krw(sw);
    }

    /*!
     * \brief The derivative of the relative permeability for the wetting phase w.r.t. saturation
     */
    Scalar dkrw_dsw(const Scalar sw) const
    {
        if (sw > swInterval_[0] && sw < swInterval_[1])
            return krwSpline_->evalDerivative(sw);

        return TwoPMaterialLaw::dkrw_dsw(sw);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     */
    Scalar krn(const Scalar sw) const
    {
        if (sw > swInterval_[0] && sw < swInterval_[1])
            return krnSpline_->eval(sw);

        return TwoPMaterialLaw::krn(sw);
    }

    /*!
     * \brief The derivative of the relative permeability for the non-wetting phase w.r.t. saturation
     */
    Scalar dkrn_dsw(const Scalar sw) const
    {
        if (sw > swInterval_[0] && sw < swInterval_[1])
            return krnSpline_->evalDerivative(sw);

        return TwoPMaterialLaw::dkrn_dsw(sw);
    }

private:
    template<class Function>
    std::unique_ptr<MonotoneCubicSpline<Scalar>>
    makeSweSpline_(const Function& f, bool invert = false) const
    {
        const auto sw = linspace(swInterval_[0], swInterval_[1], numSwSamples_);

        auto values = sw;
        std::transform(sw.begin(), sw.end(), values.begin(), f);

        if (invert)
            return std::make_unique<MonotoneCubicSpline<Scalar>>(values, sw);
        else
            return std::make_unique<MonotoneCubicSpline<Scalar>>(sw, values);
    }

    std::unique_ptr<MonotoneCubicSpline<Scalar>> pcSpline_;
    std::unique_ptr<MonotoneCubicSpline<Scalar>> krwSpline_;
    std::unique_ptr<MonotoneCubicSpline<Scalar>> krnSpline_;

    std::array<Scalar, 2> swInterval_;
    std::array<Scalar, 2> swIntervalPc_;
    std::size_t numSwSamples_;
};

} // end namespace Dumux::FluidMatrix

#endif
