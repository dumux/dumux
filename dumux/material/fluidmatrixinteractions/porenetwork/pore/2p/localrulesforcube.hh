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
 * \brief Pore-local pc-Sw curves for cubic pore bodies.
 */
#ifndef DUMUX_PNM_2P_LOCAL_RULES_FOR_CUBE_HH
#define DUMUX_PNM_2P_LOCAL_RULES_FOR_CUBE_HH

#include <cmath>
#include <dumux/porenetworkflow/common/poreproperties.hh>
#include <dumux/common/optionalscalar.hh>
#include <dumux/common/spline.hh>
#include "baselocalrules.hh"

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the simplified pore-local capillary pressure-saturation curve
 *        according to Joekar-Niasar et al., 2010.
 */
struct TwoPLocalRulesCubeJoekarNiasar
{

    /*!
     * \brief The parameter type
     * \tparam Scalar The scalar type
     */
    template<class Scalar>
    struct Params
    {
        template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
        void update(const SpatialParams& spatialParams,
                    const Element& element,
                    const SubControlVolume& scv,
                    const ElemSol& elemSol)
        {
            const auto& gridGeometry = spatialParams.gridGeometry();
            shape_ = gridGeometry.poreGeometry()[scv.dofIndex()];
            radius_ = spatialParams.poreRadius(element, scv, elemSol);

            static const Scalar surfaceTension = getParam<Scalar>("SpatialParams.SurfaceTension", 0.0725); // TODO
            surfaceTension_ = surfaceTension;
            updated_ = true;
        }

        Pore::Shape poreShape() const { return shape_; }

        Scalar poreRadius() const { return radius_; }

        Scalar surfaceTension() const { return surfaceTension_; }

        bool isUpdated() const
        { return updated_; }

        bool operator== (const Params& p) const
        {
            return Dune::FloatCmp::eq(radius_, p.radius_, 1e-6)
                   && Dune::FloatCmp::eq(surfaceTension_, p.surfaceTension_, 1e-6)
                   && shape_ == p.shape_;
        }

    private:
        Pore::Shape shape_;
        Scalar radius_;
        Scalar surfaceTension_;
        bool updated_ = false;
    };

    static constexpr bool supportsMultipleGeometries()
    { return false; }

    /*!
     * \brief The capillary pressure-saturation curve according to Joekar-Niasar et al., 2010.
     *
     *  \f$\mathrm{
     *  p_C = \frac{2*\sigma}{R(1-e^{-6.83S})}
     *  }\f$
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{S_{w,i}}\f$ at pore \f$i\f$
     * \param params The parameters container
     */
    template<class Scalar>
    static Scalar pc(const Scalar sw, const Params<Scalar>& params)
    {
        assert(0 <= sw && sw <= 1);
        assert(params.poreShape() == Pore::Shape::cube);
        const Scalar poreRadius = params.poreRadius();
        const Scalar sigma = params.surfaceTension() ;
        // TODO incorporate contact angle!!!
        using std::exp;
        return 2.0*sigma / (poreRadius*(1.0 - exp(-6.83*sw))) ;
    }

    /*!
     * \brief The wetting-phase saturation of a pore body
     *
     * \param pc The capillary pressure \f$\mathrm{p_{c,i}}\f$ at pore \f$i\f$
     * \param params The parameters container
     */
    template<class Scalar>
    static Scalar sw(const Scalar pc, const Params<Scalar>& params)
    {
        assert(params.poreShape() == Pore::Shape::cube);
        const Scalar poreRadius = params.poreRadius();
        const Scalar sigma = params.surfaceTension();
        using std::log;
        return  - 1/6.83* log(1 - 1/poreRadius * 2*sigma / pc);
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the wetting phase saturation.
     *
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{S_{w,i}}\f$ at pore \f$i\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar dpc_dsw(const Scalar sw, const Params<Scalar>& params)
    {
        assert(0 <= sw && sw <= 1);
        assert(params.poreShape() == Pore::Shape::cube);
        using std::exp;
        const Scalar sigma = params.surfaceTension();
        const Scalar poreRadius = params.poreRadius();
        const Scalar e = exp(6.83*sw);
        return -(13.66*sigma*e) / (poreRadius*(e-1.0)*(e-1.0));
    }

    /*!
     * \brief DOCU
     *
     *
     * \param pc The capillary pressure \f$\mathrm{p_{c,i}}\f$ at pore \f$i\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar dsw_dpc(const Scalar pc, const Params<Scalar>& params)
    {
        return 0; // TODO
    }
};

template<class Scalar>
class TwoPLocalRulesCubeJoekarNiasarRegularization
{
    using BaseLawParams = typename TwoPLocalRulesCubeJoekarNiasar::Params<Scalar>;
public:

    template<class S>
    struct Params
    {
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

private:
        S pcLowSw_ = 0.01;
        S pcHighSw_ = 0.99;
    };

    /*!
     * \brief The available options for regularizing the pc-SW curve at high wetting-phase saturations.
     */
    enum class HighSwRegularizationMethod
    {
        linear, spline, powerLaw
    };


    template<class MaterialLaw>
    void init(const MaterialLaw* m, const Params<Scalar>& p, const std::string& paramGroup = "")
    {
        pcLowSw_ = p.pcLowSw();
        pcHighSw_ = p.pcHighSw();

        initPcParameters_(m, pcLowSw_, pcHighSw_, paramGroup);
    }

    OptionalScalar<Scalar> pc(const Scalar sw) const
    {
        const Scalar lowSw = pcLowSw_;
        const Scalar highSw = pcHighSw_;

        // make sure that the capilary pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative is calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (sw < lowSw)
            return TwoPLocalRulesCubeJoekarNiasar::pc(lowSw, *baseLawParamsPtr_) + TwoPLocalRulesCubeJoekarNiasar::dpc_dsw(lowSw, *baseLawParamsPtr_) * (sw - lowSw);

        auto linearCurveForHighSw = [&]()
        {
            const Scalar slopeHighSw = -TwoPLocalRulesCubeJoekarNiasar::pc(highSw, *baseLawParamsPtr_) / (1.0-highSw);
            return slopeHighSw*(sw - 1.0);
        };

        if (sw <= highSw)
            return TwoPLocalRulesCubeJoekarNiasar::pc(sw, *baseLawParamsPtr_); // standard
        else if (sw <= 1.0) // regularized part below sw = 1.0
        {
            using std::pow;
            if (highSwRegularizationMethod_ == HighSwRegularizationMethod::powerLaw)
                return TwoPLocalRulesCubeJoekarNiasar::pc(highSw, *baseLawParamsPtr_) * pow(((1.0-sw)/(1.0-highSw)), 1.0/3.0);

            else if (highSwRegularizationMethod_ == HighSwRegularizationMethod::linear)
                return linearCurveForHighSw();

            else if (highSwRegularizationMethod_ == HighSwRegularizationMethod::spline)
            {
                // use spline between threshold swe and 1.0
                const Scalar yTh = TwoPLocalRulesCubeJoekarNiasar::pc(highSw, *baseLawParamsPtr_);
                // using zero derivatives at the beginning and end of the spline seems to work best
                // for some reason ...
                Spline<Scalar> sp(highSw, 1.0, // x0, x1
                                  yTh, 0, // y0, y1
                                  0.0, 0.0); // m0, m1

                return sp.eval(sw);
            }
            else
                DUNE_THROW(Dune::NotImplemented, "Regularization not method not implemented");
        }
        else // regularized part above sw = 1.0
            return linearCurveForHighSw();
    }

    OptionalScalar<Scalar> sw(const Scalar pc) const
    {
        // TODO!!!
        return {};
    }

    template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
    void updateParams(const SpatialParams& spatialParams,
                      const Element& element,
                      const SubControlVolume& scv,
                      const ElemSol& elemSol)
    {
    }


private:

    template<class MaterialLaw>
    void initPcParameters_(const MaterialLaw* m, const Scalar lowSw, const Scalar highSw, const std::string& paramGroup)
    {
        highSwRegularizationMethod_ = [&]()
        {
            const auto input = getParamFromGroup<std::string>(paramGroup, "HighSwRegularizationMethod", "Linear");
            if (input == "Linear")
                return HighSwRegularizationMethod::linear;
            else if (input == "Spline")
                return HighSwRegularizationMethod::spline;
            else if (input == "PowerLaw")
                return HighSwRegularizationMethod::powerLaw;
            else
            DUNE_THROW(Dune::InvalidStateException, input << " is not a valid regularization method");
        }();

        baseLawParamsPtr_ = &m->basicParams();
    }

    Scalar pcLowSw_, pcHighSw_;
    HighSwRegularizationMethod highSwRegularizationMethod_;
    const BaseLawParams* baseLawParamsPtr_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A default configuration for using the VanGenuchten material law
 */
template<typename Scalar = double>
using TwoPLocalRulesCubeJoekarNiasarDefault = TwoPLocalRulesBase<Scalar, TwoPLocalRulesCubeJoekarNiasar, TwoPLocalRulesCubeJoekarNiasarRegularization<Scalar>>;

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A default configuration without regularization for using the VanGenuchten material law
 */
template<typename Scalar = double>
using TwoPLocalRulesCubeJoekarNiasarNoReg = TwoPLocalRulesBase<Scalar, TwoPLocalRulesCubeJoekarNiasar, NoRegularization>;

}

#endif
