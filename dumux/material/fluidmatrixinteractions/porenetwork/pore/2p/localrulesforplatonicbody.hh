// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \ingroup PoreNetworkModels
 * \brief Pore-local pc-Sw curves for for platonic bodies
 *        (tetrahedron, cube, octahedron, dodecahedron, icosahedron).
 */
#ifndef DUMUX_PNM_2P_LOCAL_RULES_FOR_PLATONIC_BODY_HH
#define DUMUX_PNM_2P_LOCAL_RULES_FOR_PLATONIC_BODY_HH

#include <algorithm>
#include <cassert>
#include <cmath>
#include <optional>
#include <dune/common/math.hh>
#include <dumux/common/optionalscalar.hh>
#include <dumux/common/spline.hh>
#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/porenetwork/common/poreproperties.hh>
#include "singleshapelocalrules.hh"

namespace Dumux::PoreNetwork::FluidMatrix {

/*!
 * \brief The parameter type
 * \tparam Scalar The scalar type
 */
template<class Scalar>
struct PlatonicBodyParams
{
    PlatonicBodyParams() = default;

    PlatonicBodyParams& setPoreInscribedRadius(Scalar r) { radius_ = r; return *this;}
    PlatonicBodyParams& setPoreShape(Pore::Shape s) { shape_ = s; return *this;}
    PlatonicBodyParams& setSurfaceTension(Scalar st) { surfaceTension_ = st; return *this;}

    template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
    PlatonicBodyParams(const SpatialParams& spatialParams,
                       const Element& element,
                       const SubControlVolume& scv,
                       const ElemSol& elemSol)
    : shape_(spatialParams.gridGeometry().poreGeometry()[scv.dofIndex()])
    , radius_(spatialParams.poreInscribedRadius(element, scv, elemSol))
    , surfaceTension_(spatialParams.surfaceTension(element, scv, elemSol))
    {}

    template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
    void update(const SpatialParams& spatialParams,
                const Element& element,
                const SubControlVolume& scv,
                const ElemSol& elemSol)
    {
        const auto& gridGeometry = spatialParams.gridGeometry();
        shape_ = gridGeometry.poreGeometry()[scv.dofIndex()];
        radius_ = spatialParams.poreInscribedRadius(element, scv, elemSol);
        surfaceTension_ = spatialParams.surfaceTension(element, scv, elemSol);
    }

    Pore::Shape poreShape() const { return shape_; }

    Scalar poreInscribedRadius() const { return radius_; }

    Scalar surfaceTension() const { return surfaceTension_; }

    bool operator== (const PlatonicBodyParams& p) const
    {
        return Dune::FloatCmp::eq(radius_, p.radius_, 1e-6)
               && Dune::FloatCmp::eq(surfaceTension_, p.surfaceTension_, 1e-6)
               && shape_ == p.shape_;
    }

private:
    Pore::Shape shape_;
    Scalar radius_;
    Scalar surfaceTension_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \ingroup PoreNetworkModels
 * \brief Implementation of the simplified pore-local capillary pressure-saturation curve
 *        for platonic bodies (tetrahedron, cube, octahedron, dodecahedron, icosahedron).
 *
 *        See Joekar-Niasar et al., 2010 and Sweijen et al., 2018.
 */
template<Pore::Shape shape>
struct TwoPLocalRulesPlatonicBody
{
    template<class Scalar>
    using Params = PlatonicBodyParams<Scalar>;

    template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
    static auto makeParams(const SpatialParams& spatialParams,
                           const Element& element,
                           const SubControlVolume& scv,
                           const ElemSol& elemSol)
    {
        using Scalar = std::decay_t<decltype(spatialParams.poreInscribedRadius(element, scv, elemSol))>;
        return Params<Scalar>(spatialParams, element, scv, elemSol);
    }

    /*!
     * \brief The capillary pressure-saturation curve
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{S_{w,i}}\f$ at pore \f$i\f$
     * \param params The parameters container
     */
    template<class Scalar>
    static Scalar pc(Scalar sw, const Params<Scalar>& params)
    {
        assert(isPlatonicBody_(params.poreShape()));

        using std::clamp;
        sw = clamp(sw, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar poreRadius = params.poreInscribedRadius();
        const Scalar sigma = params.surfaceTension();
        // TODO incorporate contact angle!!!

        using std::exp;
        return 2.0*sigma / (poreRadius*(1.0 - exp(-expFactor_<Scalar>()*sw)));
    }

    /*!
     * \brief The wetting-phase saturation of a pore body
     *
     * \param pc The capillary pressure \f$\mathrm{p_{c,i}}\f$ at pore \f$i\f$
     * \param params The parameters container
     */
    template<class Scalar>
    static Scalar sw(Scalar pc, const Params<Scalar>& params)
    {
        assert(isPlatonicBody_(params.poreShape()));

        using std::max;
        pc = max(pc, 0.0); // the equation below is undefined for negative pcs

        const Scalar poreRadius = params.poreInscribedRadius();
        const Scalar sigma = params.surfaceTension();

        using std::log;
        return -1.0/expFactor_<Scalar>()* log(1.0 - 2.0*sigma/(poreRadius*pc));
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the wetting phase saturation.
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{S_{w,i}}\f$ at pore \f$i\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar dpc_dsw(Scalar sw, const Params<Scalar>& params)
    {
        assert(isPlatonicBody_(params.poreShape()));

        using std::clamp;
        sw = clamp(sw, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar sigma = params.surfaceTension();
        const Scalar poreRadius = params.poreInscribedRadius();
        using std::exp;
        const Scalar e = exp(expFactor_<Scalar>()*sw);
        return -(2.0*expFactor_<Scalar>()*sigma*e) / (poreRadius*(1.0-e)*(1.0-e));
    }

    /*!
     * \brief The partial derivative of the wetting phase saturation
     *        w.r.t. the capillary pressure.
     *
     * \param pc The capillary pressure \f$\mathrm{p_{c,i}}\f$ at pore \f$i\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar dsw_dpc(Scalar pc, const Params<Scalar>& params)
    {
        assert(isPlatonicBody_(params.poreShape()));

        using std::max;
        pc = max(pc, 0.0); // the equation below is undefined for negative pcs

        const Scalar sigma = params.surfaceTension();
        const Scalar poreRadius = params.poreInscribedRadius();
        return sigma / (expFactor_<Scalar>()*sigma*pc - 0.5*expFactor_<Scalar>()*poreRadius*pc*pc);
    }

private:

    template<class Scalar>
    static constexpr Scalar expFactor_()
    {
        if constexpr (shape == Pore::Shape::tetrahedron)
            return 3.87;
        else if constexpr (shape == Pore::Shape::cube)
            return 6.83;
        else if constexpr (shape == Pore::Shape::octahedron)
            return 8.71;
        else if constexpr (shape == Pore::Shape::dodecahedron)
            return 22.87;
        else if constexpr (shape == Pore::Shape::icosahedron)
            return 24.11;
        else
        {
            static_assert(AlwaysFalse<Scalar>::value, "Shape not supported");
            return 0;
        }
    }

    static constexpr bool isPlatonicBody_(Pore::Shape s)
    {
        return s == Pore::Shape::tetrahedron ||
               s == Pore::Shape::cube ||
               s == Pore::Shape::octahedron ||
               s == Pore::Shape::dodecahedron ||
               s == Pore::Shape::icosahedron;
    }
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \ingroup PoreNetworkModels
 * \brief Two-phase rules for regularizing the pc-SW for platonic bodies.
 */
template<class Scalar, class BaseLaw>
class TwoPLocalRulesPlatonicBodyRegularization
{
    using ThisType = TwoPLocalRulesPlatonicBodyRegularization<Scalar, BaseLaw>;
    using BaseLawParams = typename BaseLaw::template Params<Scalar>;
public:

    /*!
     * \brief The available options for regularizing the pc-SW curve at high wetting-phase saturations.
     */
    enum class HighSwRegularizationMethod
    {
        linear, spline, powerLaw
    };

    template<class S>
    struct Params
    {
        // export method type
        using HighSwRegularizationMethod = ThisType::HighSwRegularizationMethod;
        /*!
         * \brief Set the threshold saturation below which the capillary pressure is regularized.
         *
         * Most problems are very sensitive to this value (e.g. making it smaller might
         * result in very high capillary pressures)
         */
        void setpcLowSw(S pcLowSw)
        { pcLowSw_ = pcLowSw; }

        /*!
         * \brief Threshold saturation below which the capillary pressure is regularized.
         */
        S pcLowSw() const
        { return pcLowSw_; }

        /*!
         * \brief Set the threshold saturation above which the capillary pressure is regularized.
         */
        void setpcHighSw(S pcHighSw)
        { pcHighSw_ = pcHighSw; }

        /*!
         * \brief Threshold saturation above which the capillary pressure is regularized.
         *
         * Most problems are very sensitive to this value (e.g. making it smaller might
         * result in negative capillary pressures).
         */
        S pcHighSw() const
        { return pcHighSw_; }

        /*!
         * \brief Set the regularization method for high saturations.
         */
        void setHighSwRegularizationMethod(HighSwRegularizationMethod method)
        { highSwRegularizationMethod_ = method; }

        /*!
         * \brief Return the regularization method for high saturations.
         */
        HighSwRegularizationMethod highSwRegularizationMethod() const
        { return highSwRegularizationMethod_; }

private:
        S pcLowSw_ = 0.01;
        S pcHighSw_ = 0.99;
        HighSwRegularizationMethod highSwRegularizationMethod_ = HighSwRegularizationMethod::linear;
    };


    template<class MaterialLaw>
    void init(const MaterialLaw* m, const Params<Scalar>& p, const std::string& paramGroup = "")
    {
        initPcParameters_(m, p, paramGroup);
    }

    OptionalScalar<Scalar> pc(const Scalar sw) const
    {
        // make sure that the capilary pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative is calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (sw < pcLowSw_)
            return pcLowSwPcValue_() + pcDerivativeLowSw_() * (sw - pcLowSw_);

        if (sw <= pcHighSw_)
            return {}; // standard
        else if (sw < 1.0) // regularized part below sw = 1.0
        {
            using std::pow;
            if (highSwRegularizationMethod_ == HighSwRegularizationMethod::powerLaw)
                return pcHighSwPcValue_() * pow(((1.0-sw)/(1.0-pcHighSw_)), 1.0/3.0);

            else if (highSwRegularizationMethod_ == HighSwRegularizationMethod::linear)
            {
                const Scalar slope = -pcHighSwPcValue_() / (1.0 - pcHighSw_);
                return pcHighSwPcValue_() + (sw - pcHighSw_) * slope;
            }

            else if (highSwRegularizationMethod_ == HighSwRegularizationMethod::spline)
                return pcSpline_().eval(sw);

            else
                DUNE_THROW(Dune::NotImplemented, "Regularization not method not implemented");
        }
        else // regularized part above sw = 1.0
            return pcDerivativeHighSwEnd_()*(sw - 1.0);
    }

    /*!
     * \brief The regularized saturation-capillary pressure curve
     */
    OptionalScalar<Scalar> sw(const Scalar pc) const
    {
        if (pc <= 0.0)
        {
            if (pcHighSw_ >= 1.0)
                return 1.0; // no regularization for high sw
            else
                return pc/pcDerivativeHighSwEnd_() + 1.0;
        }

        // low saturation
        if (pc > pcLowSwPcValue_())
            return (pc - pcLowSwPcValue_())/pcDerivativeLowSw_() + pcLowSw_;

        // high saturation
        else if (pc <= pcHighSwPcValue_())
        {
            if (highSwRegularizationMethod_ == HighSwRegularizationMethod::powerLaw)
            {
                // invert power law
                using Dune::power;
                return power(pc/pcHighSwPcValue_(), 3) * (pcHighSw_ - 1.0) + 1.0;
            }

            else if (highSwRegularizationMethod_ == HighSwRegularizationMethod::linear)
                return pc/pcDerivativeHighSwEnd_() + 1.0;

            else if (highSwRegularizationMethod_ == HighSwRegularizationMethod::spline)
                // invert spline
                return pcSpline_().intersectInterval(pcHighSw_, 1.0, 0.0, 0.0, 0.0, pc);
        }

        // no regularization
        return {};
    }

    /*!
     * \brief The regularized partial derivative of the capillary pressure w.r.t. the saturation
     */
    OptionalScalar<Scalar> dpc_dsw(const Scalar sw) const
    {
        if (sw <= pcLowSw_)
            return pcDerivativeLowSw_();

        else if (sw >= 1.0)
            return pcDerivativeHighSwEnd_();

        else if (sw > pcHighSw_)
        {
            if (highSwRegularizationMethod_ == HighSwRegularizationMethod::powerLaw)
            {
                using std::pow;
                return pcHighSwPcValue_()/3.0 * 1.0 /((pcHighSw_-1.0) * pow((sw-1.0)/(pcHighSw_-1.0), 2.0/3.0));
            }

            else if (highSwRegularizationMethod_ == HighSwRegularizationMethod::linear)
                return pcDerivativeHighSwEnd_();

            else
                return pcSpline_().evalDerivative(sw);
        }

        else
            return {}; // no regularization
    }

    /*!
     * \brief The regularized partial derivative of the saturation to the capillary pressure
     */
    OptionalScalar<Scalar> dsw_dpc(const Scalar pc) const
    {
        if (pc <= 0.0)
        {
            if (pcHighSw_ >= 1.0)
                return 0.0;
            else
                return 1.0/pcDerivativeHighSwEnd_();
        }

        // derivative of the inverse of the function is one over derivative of the function
        else if (pc <= pcHighSwPcValue_())
        {
            if (highSwRegularizationMethod_ == HighSwRegularizationMethod::powerLaw)
            {
                using Dune::power;
                return (3.0*pcHighSw_ - 3.0) * power(pc, 2) / power(pcHighSwPcValue_(), 3);
            }

            else if (highSwRegularizationMethod_ == HighSwRegularizationMethod::linear)
                return 1.0/pcDerivativeHighSwEnd_();

            else
                return 1.0/pcSpline_().evalDerivative(pcSpline_().intersectInterval(pcHighSw_, 1.0, 0.0, 0.0, 0.0, pc));
        }

        else if (pc >= pcLowSwPcValue_())
            return 1.0/pcDerivativeLowSw_();

        else
            return {}; // no regularization
    }

    template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
    void updateParams(const SpatialParams& spatialParams,
                      const Element& element,
                      const SubControlVolume& scv,
                      const ElemSol& elemSol)
    {}

private:

    template<class MaterialLaw>
    void initPcParameters_(const MaterialLaw* m, const Params<Scalar>& params, const std::string& paramGroup)
    {
        highSwRegularizationMethod_ = params.highSwRegularizationMethod();

        // maybe overwrite method using input
        static const bool hasMethod = hasParamInGroup(paramGroup, "HighSwRegularizationMethod");
        if (hasMethod)
        {
            // set method only once
            static const HighSwRegularizationMethod inputmethod = [&]
            {
                static const auto input = getParamFromGroup<std::string>(paramGroup, "HighSwRegularizationMethod");
                if (input == "Linear")
                    return HighSwRegularizationMethod::linear;
                else if (input == "Spline")
                    return HighSwRegularizationMethod::spline;
                else if (input == "PowerLaw")
                    return HighSwRegularizationMethod::powerLaw;
                else
                    DUNE_THROW(Dune::InvalidStateException, input << " is not a valid regularization method");
            }();

            highSwRegularizationMethod_ = inputmethod;
        }

        static const bool splineZeroSlope = getParamFromGroup<bool>(paramGroup, "HighSwSplineZeroSlope", true);
        highSwSplineZeroSlope_ = splineZeroSlope;

        // print name only once
        [[maybe_unused]] static const bool printName = [&]
        {
            std::string name;
            if (highSwRegularizationMethod_ == HighSwRegularizationMethod::linear)
                name = "Linear";
            else if (highSwRegularizationMethod_ == HighSwRegularizationMethod::powerLaw)
                name =  "PowerLaw";
            else if (highSwRegularizationMethod_ == HighSwRegularizationMethod::spline)
                name = "Spline";
            else
                DUNE_THROW(Dune::NotImplemented, "Regularization not method not implemented");

            std::cout << "\n*****\nUsing " << name << " as regularization method for high Sw\n*****" << std::endl;
            return true;
        }();

        using std::isnan;
        static const auto pcLowSwInput = getParamFromGroup<Scalar>(paramGroup, "RegularizationLowSw", params.pcLowSw());
        pcLowSw_ = !isnan(pcLowSwInput) ? pcLowSwInput : params.pcLowSw();

        static const auto pcHighSwInput = getParamFromGroup<Scalar>(paramGroup, "RegularizationHighSw", std::numeric_limits<Scalar>::quiet_NaN());
        pcHighSw_ = !isnan(pcHighSwInput) ? pcHighSwInput : params.pcHighSw();

        baseLawParamsPtr_ = &m->basicParams();

        static const auto highSwFixedSlopeInput = getParamFromGroup<Scalar>(paramGroup, "RegularizationHighSwFixedSlope", std::numeric_limits<Scalar>::quiet_NaN());
        highSwFixedSlope_ = highSwFixedSlopeInput;

        // Note: we do not pre-calculate all end-point values (e.g.,  pcLowSwPcValue_ and pcDerivativeLowSw_)
        // as done in, e.g., VanGenuchten. This is because this law will generally be instantiated only as
        // a temporary object since all pore bodies generally have different parameters.
        // We still cache above-mentioned values, but only if they are actually needed (see below).
    }

    // the capillary pressure for the lower saturation threshold
    Scalar pcLowSwPcValue_() const
    {
        // calculated value within first function call, used cached value otherwise
        if (!optionalPcLowSwPcValue_)
            optionalPcLowSwPcValue_ = BaseLaw::pc(pcLowSw_, *baseLawParamsPtr_);
        return optionalPcLowSwPcValue_.value();
    }

    // dpc_dsw for the lower saturation threshold
    Scalar pcDerivativeLowSw_() const
    {
        if (!optionalPcDerivativeLowSw_)
            optionalPcDerivativeLowSw_ = BaseLaw::dpc_dsw(pcLowSw_, *baseLawParamsPtr_);
        return optionalPcDerivativeLowSw_.value();
    }

    // the capillary pressure for the upper saturation threshold
    Scalar pcHighSwPcValue_() const
    {
        if (!optionalPcHighSwPcValue_)
            optionalPcHighSwPcValue_ = BaseLaw::pc(pcHighSw_, *baseLawParamsPtr_);
        return optionalPcHighSwPcValue_.value();
    }

    // dpc_dsw at Sw = 1.0
    Scalar pcDerivativeHighSwEnd_() const
    {
        using std::isnan;
        if (!isnan(highSwFixedSlope_))
            return highSwFixedSlope_;
        else
            return (0.0 - pcHighSwPcValue_()) / (1.0 - pcHighSw_);
    }

    const Spline<Scalar>& pcSpline_() const
    {
        if (!optionalPcSpline_)
        {
            const auto slopes = highSwSplineZeroSlope_ ? std::array{0.0, 0.0}
                                                       : std::array{BaseLaw::dpc_dsw(pcHighSw_, *baseLawParamsPtr_), pcDerivativeHighSwEnd_()};

            optionalPcSpline_ = Spline<Scalar>(pcHighSw_, 1.0, // x0, x1
                                               pcHighSwPcValue_(), 0, // y0, y1
                                               slopes[0], slopes[1]); // m0, m1
        }
        return optionalPcSpline_.value();
    }

    Scalar pcLowSw_, pcHighSw_;
    mutable OptionalScalar<Scalar> optionalPcLowSwPcValue_, optionalPcDerivativeLowSw_;
    mutable OptionalScalar<Scalar> optionalPcHighSwPcValue_;
    HighSwRegularizationMethod highSwRegularizationMethod_;
    const BaseLawParams* baseLawParamsPtr_;
    mutable std::optional<Spline<Scalar>> optionalPcSpline_;
    bool highSwSplineZeroSlope_;
    Scalar highSwFixedSlope_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A default configuration for using the VanGenuchten material law
 */
template<Pore::Shape shape, typename Scalar = double>
using TwoPLocalRulesPlatonicBodyDefault = SingleShapeTwoPLocalRules<Scalar,
                                                                    TwoPLocalRulesPlatonicBody<shape>,
                                                                    TwoPLocalRulesPlatonicBodyRegularization<Scalar,
                                                                                                             TwoPLocalRulesPlatonicBody<shape>>>;

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A default configuration without regularization for using the VanGenuchten material law
 */
template<Pore::Shape shape, typename Scalar = double>
using TwoPLocalRulesPlatonicBodyNoReg = SingleShapeTwoPLocalRules<Scalar, TwoPLocalRulesPlatonicBody<shape>, Dumux::FluidMatrix::NoRegularization>;

} // end namespace Dumux::FluidMatrix

#endif
