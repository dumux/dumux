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
 * \brief Implementation of the capillary pressure / relPerm <-> saturation relation
 *        using a linear relation smoothed at the upper and lower bounds for kr.
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWOP_SMOOTHED_LINEAR_LAW_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_TWOP_SMOOTHED_LINEAR_LAW_HH


#include <cmath>
#include <algorithm>

#include <dumux/common/parameters.hh>
#include <dumux/common/spline.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabsdefaultpolicy.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implements a linear saturation-capillary pressure relation
 *
 * The entry pressure is reached at \f$\mathrm{\overline{S}_w = 1}\f$, the maximum
 * capillary pressure is observed at \f$\mathrm{\overline{S}_w = 0}\f$.
 *
 * The relative permeabilities are 0 or 1 outside of the range of effective saturation.
 * However, the transition between the linearly changing and the constant part is not smooth but with a kink.
 * The Newton scheme does not like that. Therefore a smooth transition is accomplished by interpolating these
 * regions with a spline.
 *
 * An example of the regularization of the relative permeability is shown below:
 * \image html regularizedLinearKr.png
 */
template<class ScalarType, class EffToAbsPolicy = TwoPEffToAbsDefaultPolicy>
class SmoothedLinearLaw : public Adapter<SmoothedLinearLaw<ScalarType, EffToAbsPolicy>, PcKrSw>
{

public:
    using Scalar = ScalarType;
    using EffToAbsParams = typename EffToAbsPolicy::template Params<Scalar>;
    using EffToAbs = EffToAbsPolicy;

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
        Params(Scalar pe, Scalar pcMax, Scalar krLowS, Scalar krHighS)
        : pe_(pe), pcMax_(pcMax), krLowS_(krLowS), krHighS_(krHighS)
        {}

        Scalar pe() const { return pe_; }
        void setPe(Scalar pe) { pe_ = pe; }

        Scalar pcMax() const { return pcMax_; }
        void setPcMax(Scalar pcMax) { pcMax_ = pcMax; }

        Scalar krLowS() const { return krLowS_; }
        void setKrLowS(Scalar krLowS) { krLowS_ = krLowS; }

        Scalar krHighS() const { return krHighS_; }
        void setKrHighS(Scalar krHighS) { krHighS_ = krHighS; }

        bool operator== (const Params& p) const
        {
            return Dune::FloatCmp::eq(pe(), p.pe(), 1e-6)
                   && Dune::FloatCmp::eq(pcMax(), p.pcMax(), 1e-6)
                   && Dune::FloatCmp::eq(krLowS(), p.krLowS(), 1e-6)
                   && Dune::FloatCmp::eq(krHighS(), p.krHighS(), 1e-6);
        }

    private:
        Scalar pe_, pcMax_, krLowS_, krHighS_;
    };

    /*!
     * \brief Deleted default constructor (so we are never in an undefined state)
     * \note store owning pointers to laws instead if you need default-constructible objects
     */
    SmoothedLinearLaw() = delete;

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    explicit SmoothedLinearLaw(const std::string& paramGroup)
    : SmoothedLinearLaw(makeParams(paramGroup), EffToAbs::template makeParams<Scalar>(paramGroup))
    {}

    /*!
     * \brief Construct from parameter structs
     * \note More efficient constructor but you need to ensure all parameters are initialized
     */
    SmoothedLinearLaw(const Params& params,
                      const EffToAbsParams& effToAbsParams = {})
    : params_(params)
    , effToAbsParams_(effToAbsParams)
    , splineM_((1.0 - ((1.0 - params_.krHighS()) + params_.krLowS())/2.0 )
               / (1.0 - (1.0 - params_.krHighS()) - params_.krLowS()))
    , splineLowS_(0.0, params_.krLowS(), // x1, x2
                  0.0, params_.krLowS()/2.0, // y1, y2
                  0.0, splineM_) // m1, m2
    , splineHighS_(params_.krHighS(), 1.0, // x1, x2
                   1.0 - (1.0 - params_.krHighS())/2.0, 1.0, // y1, y2
                   splineM_, 0.0) // m1, m2
    {}

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    static Params makeParams(const std::string& paramGroup)
    {
        const auto pe = getParamFromGroup<Scalar>(paramGroup, "SmoothedLinearLawPe");
        const auto pcMax = getParamFromGroup<Scalar>(paramGroup, "SmoothedLinearLawPcMax");
        const auto krLowS = getParamFromGroup<Scalar>(paramGroup, "SmoothedLinearLawKrLowS");
        const auto krHighS = getParamFromGroup<Scalar>(paramGroup, "SmoothedLinearLawKrHighS");
        return {pe, pcMax, krLowS, krHighS};
    }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    Scalar pc(Scalar swe) const
    {
        return (1.0 - swe)*(params_.pcMax() - params_.pe()) + params_.pe();
    }

    /*!
     * \brief The inverse saturation-capillary pressure curve
     */
    Scalar swe(Scalar pc) const
    {
        return 1.0 - (pc - params_.pe())/(params_.pcMax() - params_.pe());
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     */
    Scalar endPointPc() const
    { return params_.pe(); }

    /*!
     * \brief The partial derivative of the capillary pressure w.r.t. the effective saturation
     */
    Scalar dpc_dswe(Scalar swe) const
    {
        return params_.pe() - params_.pcMax();
    }

    /*!
     * \brief The partial derivative of the effective saturation w.r.t. the capillary pressure
     */
    Scalar dswe_dpc(Scalar pc) const
    {
        return 1.0/(params_.pe() - params_.pcMax());
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     */
    Scalar krw(Scalar swe) const
    {
        return relperm_(swe);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase.
     */
    Scalar krn(Scalar swe) const
    {
        Scalar sne = 1.0 - swe;
        return relperm_(sne);
    }

    /*!
     * \brief Equality comparison with another instance
     */
    bool operator== (const SmoothedLinearLaw<Scalar, EffToAbs>& o) const
    {
        return params_ == o.params_
               && effToAbsParams_ == o.effToAbsParams_;
    }

    const EffToAbsParams& effToAbsParams() const
    { return effToAbsParams_; }

private:

    Scalar relperm_(Scalar S) const
    {
        const Scalar lowS = params_.krLowS();
        const Scalar highS = params_.krHighS();

        // check whether the saturation is unpyhsical
        if (S >= 1.0)
            return 1.0;
        else if (S <= 0.0)
            return 0;
        // check wether the permeability needs to be regularized
        else if (S < lowS)
            return splineLowS_.eval(S);
        else if (S > highS)
            return splineHighS_.eval(S);

        // straight line for S \in [lowS, highS]
        return lowS/2.0 + splineM_*(S - lowS);
    }

    Params params_;
    EffToAbsParams effToAbsParams_;
    Scalar splineM_;
    Dumux::Spline<Scalar> splineLowS_;
    Dumux::Spline<Scalar> splineHighS_;
};
} // end namespace Dumux::FluidMatrix

#endif
