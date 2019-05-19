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
 * \brief   Implementation of the capillary pressure and
 *          relative permeability <-> saturation relations for two-phase models
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWOP_MATERIAL_LAW_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_TWOP_MATERIAL_LAW_HH

#include <dumux/common/parameters.hh>
#include <dumux/material/fluidmatrixinteractions/2pnew/efftoabsdefaultpolicy.hh>
#include <dumux/material/fluidmatrixinteractions/2pnew/regularization.hh>

namespace Dumux {
namespace FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief   Implementation of the capillary pressure and
 *          relative permeability <-> saturation relations according to van Genuchten.
 * \tparam ScalarType the scalar type
 * \tparam BaseLaw the base law (e.g. VanGenuchten, BrooksCorey, Linear, ...)
 * \tparam Regularization the regularization type (set to NoRegularization to turn it off)
 * \tparam EffToAbsPolicy the policy how to convert effective <-> absolute saturations
 */
template<class ScalarType,
         class BaseLaw,
         class Regularization = TwoPDefaultRegularization<ScalarType>,
         class EffToAbsPolicy = TwoPEffToAbsDefaultPolicy>
class TwoPMaterialLaw
{
    using NoRegularization = NoTwoPRegularization<ScalarType>;
public:

    using Scalar = ScalarType;
    using BaseLawParams = typename BaseLaw::template Params<Scalar>;
    using EffToAbsParams = typename EffToAbsPolicy::template Params<Scalar>;
    using RegularizationParams = typename Regularization::template Params<Scalar>;

    /*!
     * \brief Return whether this law is regularized
     */
    static constexpr bool isRegularized()
    { return !std::is_same<Regularization, NoRegularization>::value; }

    /*!
     * \brief Deleted default constructor (so we are never in an undefined state)
     * \note store pointers to laws instead
     */
    TwoPMaterialLaw() = delete;

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    explicit TwoPMaterialLaw(const std::string& paramGroup)
    {
        baseParams_ = BaseLaw::template makeParams<Scalar>(paramGroup);
        effToAbsParams_ = EffToAbsPolicy::template makeParams<Scalar>(paramGroup);
        regularization_.init(this, paramGroup);
    }

    /*!
     * \brief Construct from parameter structs
     * \note More efficient constructor but you need to ensure all parameters are initialized
     */
    TwoPMaterialLaw(const BaseLawParams& baseParams,
                    const EffToAbsParams& effToAbsParams = {},
                    const RegularizationParams& regParams = {})
    : baseParams_(baseParams)
    , effToAbsParams_(effToAbsParams)
    {
        regularization_.init(this, baseParams, effToAbsParams, regParams);
    }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    template<bool enableRegularization = isRegularized()>
    Scalar pc(const Scalar sw) const
    {
        const auto swe = EffToAbsPolicy::swToSwe(sw, effToAbsParams_);
        if (enableRegularization)
        {
            const auto regularized = regularization_.pc(swe);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::pc(swe, baseParams_);
    }

    /*!
     * \brief The partial derivative of the capillary pressure w.r.t. the saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar dpc_dsw(const Scalar sw) const
    {
        const auto swe = EffToAbsPolicy::swToSwe(sw, effToAbsParams_);
        if (enableRegularization)
        {
            const auto regularized = regularization_.dpc_dsw(swe);
            if (regularized)
                return regularized.value()*EffToAbsPolicy::dswe_dsw(effToAbsParams_);
        }

        return BaseLaw::dpc_dswe(swe, baseParams_)*EffToAbsPolicy::dswe_dsw(effToAbsParams_);
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     */
    Scalar endPointPc() const
    {
        return BaseLaw::endPointPc(baseParams_);
    }

    /*!
     * \brief The saturation-capillary pressure curve
     */
    template<bool enableRegularization = isRegularized()>
    Scalar sw(const Scalar pc) const
    {
        if (enableRegularization)
        {
            const auto regularized = regularization_.sw(pc);
            if (regularized)
                return EffToAbsPolicy::sweToSw(regularized.value(), effToAbsParams_);
        }

        return EffToAbsPolicy::sweToSw(BaseLaw::sw(pc, baseParams_), effToAbsParams_);
    }

    /*!
     * \brief The partial derivative of the saturation to the capillary pressure
     */
    template<bool enableRegularization = isRegularized()>
    Scalar dsw_dpc(const Scalar pc) const
    {
        if (enableRegularization)
        {
            const auto regularized = regularization_.dsw_dpc(pc);
            if (regularized)
                return regularized.value()*EffToAbsPolicy::dsw_dswe(effToAbsParams_);
        }

        return BaseLaw::dswe_dpc(pc, baseParams_)*EffToAbsPolicy::dsw_dswe(effToAbsParams_);
    }

    /*!
     * \brief The relative permeability for the wetting phase
     */
    template<bool enableRegularization = isRegularized()>
    Scalar krw(const Scalar sw) const
    {
        const auto swe = EffToAbsPolicy::swToSwe(sw, effToAbsParams_);
        if (enableRegularization)
        {
            const auto regularized = regularization_.krw(swe);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::krw(swe, baseParams_);
    }

    /*!
     * \brief The derivative of the relative permeability for the wetting phase w.r.t. saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar dkrw_dsw(const Scalar sw) const
    {
        const auto swe = EffToAbsPolicy::swToSwe(sw, effToAbsParams_);
        if (enableRegularization)
        {
            const auto regularized = regularization_.dkrw_dswe(swe);
            if (regularized)
                return regularized.value()*EffToAbsPolicy::dswe_dsw(effToAbsParams_);
        }

        return BaseLaw::dkrw_dswe(swe, baseParams_)*EffToAbsPolicy::dswe_dsw(effToAbsParams_);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     */
    template<bool enableRegularization = isRegularized()>
    Scalar krn(const Scalar sw) const
    {
        const auto swe = EffToAbsPolicy::swToSwe(sw, effToAbsParams_);
        if (enableRegularization)
        {
            const auto regularized = regularization_.krn(swe);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::krw(swe, baseParams_);
    }

    /*!
     * \brief The derivative of the relative permeability for the non-wetting phase w.r.t. saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar dkrn_dsw(const Scalar sw) const
    {
        const auto swe = EffToAbsPolicy::swToSwe(sw, effToAbsParams_);
        if (enableRegularization)
        {
            const auto regularized = regularization_.dkrn_dswe(swe);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::dkrn_dswe(swe, baseParams_)*EffToAbsPolicy::dswe_dsw(effToAbsParams_);
    }

private:
    BaseLawParams baseParams_;
    EffToAbsParams effToAbsParams_;
    Regularization regularization_;
};

} // end namespace FluidMatrix
} // end namespace Dumux

#endif
