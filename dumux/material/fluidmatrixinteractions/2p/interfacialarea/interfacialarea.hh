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
 * \brief   Implementation helper for capillary-pressure-saturation-interfacial-area relations
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWO_P_INTERFACIAL_AREA_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_TWO_P_INTERFACIAL_AREA_HH

#include <dumux/common/parameters.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabsdefaultpolicy.hh>
#include <dumux/material/fluidmatrixinteractions/2p/noregularization.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Wrapper class to implement regularized laws (pc-sw-a)
 *        with a conversion policy between absolution and effective saturations
 * \tparam ScalarType the scalar type
 * \tparam BaseLaw the base law (e.g. VanGenuchten, BrooksCorey, Linear, ...)
 * \tparam Regularization the regularization type (set to NoAwnRegularization to turn it off)
 * \tparam EffToAbsPolicy the policy how to convert effective <-> absolute saturations
 *
 * \note The regularization interface is expected to return Dumux::OptionalScalars which
 *       are wrappers around a Scalar type that provide a boolean operator to
 *       check whether the result is valid. If the regularization returns
 *       a non-valid value, it means that the given parameter
 *       range is outside the regularized region.
 *       For that case we forward to the call to the standard law.
 */
template<class ScalarType,
         class BaseLaw,
         template<class> class InterfaceType,
         class Regularization = NoRegularization,
         class EffToAbsPolicy = TwoPEffToAbsDefaultPolicy>
class InterfacialArea : public Adapter<InterfacialArea<ScalarType, BaseLaw, InterfaceType, Regularization, EffToAbsPolicy>, InterfaceType>
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
    InterfacialArea() = delete;

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    explicit InterfacialArea(const std::string& paramGroup)
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
    InterfacialArea(const BasicParams& baseParams,
                    const EffToAbsParams& effToAbsParams = {},
                    const RegularizationParams& regParams = {})
    : basicParams_(baseParams)
    , effToAbsParams_(effToAbsParams)
    {
        if constexpr (isRegularized())
            regularization_.init(this, baseParams, effToAbsParams, regParams);
    }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    template<bool enableRegularization = isRegularized()>
    Scalar area(const Scalar sw, const Scalar pc) const
    {
        const auto swe = EffToAbs::swToSwe(sw, effToAbsParams_);
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.area(swe, pc);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::area(swe, pc, basicParams_);
    }

    /*!
     * \brief The partial derivative of the capillary pressure w.r.t. the saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar darea_dpc(const Scalar sw, const Scalar pc) const
    {
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.darea_dpc(sw, pc);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::darea_dpc(sw, pc, basicParams_);
    }

    /*!
     * \brief The partial derivative of the saturation to the capillary pressure
     */
    template<bool enableRegularization = isRegularized()>
    Scalar darea_dsw(const Scalar sw, const Scalar pc) const
    {
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.darea_dsw(sw, pc);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::darea_dsw(sw, pc, basicParams_);
    }

    /*!
     * \brief Equality comparison with another instance
     */
    bool operator== (const InterfacialArea& o) const
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

} // end namespace Dumux::FluidMatrix

#endif
