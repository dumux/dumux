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
 *
 * \brief Implementation of a regularized version of the pore network
 *        capillary pressure / relative permeability  <-> saturation relation.
 */
#ifndef REGULARIZED_PNM_2P_LOCAL_RULES_FOR_OCTAHEDRON_HH
#define REGULARIZED_PNM_2P_LOCAL_RULES_FOR_OCTAHEDRON_HH

#include "../baselocalrules.hh"
#include "../localrules/localrulesforoctahedron.hh"

#include <dumux/common/parameters.hh>
#include <dumux/common/spline.hh>

namespace Dumux
{

/*!\ingroup fluidmatrixinteractionslaws
 *
 * \brief Implementation of the regularized  pore network
 *        capillary pressure / relative permeability  <-> saturation relation.
 *        This class bundles the "raw" curves as
 *        static members and doesn't concern itself converting
 *        absolute to effective saturations and vice versa.
 *
 *        In order to avoid very steep gradients the marginal values are "regularized".
 *        This means that in stead of following the curve of the material law in these regions, some linear approximation is used.
 *        Doing this is not worse than following the material law. E.g. for very low wetting phase values the material
 *        laws predict infinite values for \f$\mathrm{p_c}\f$ which is completely unphysical. In case of very high wetting phase
 *        saturations the difference between regularized and "pure" material law is not big.
 *
 *        Regularizing has the additional benefit of being numerically friendly: Newton's method does not like infinite gradients.
 *
 *        The implementation is accomplished as follows:
 *        - check whether we are in the range of regularization
 *         - yes: use the regularization
 *         - no: forward to the standard material law.
 *
 *         For an example figure of the regularization: RegularizedVanGenuchten
 *
 * \see PNMLocalRules
 */
template<class ScalarT>
class RegularizedTwoPLocalRulesOctahedron : public RegularizedTwoPLocalRulesBase
{
    using LocalRules = TwoPLocalRulesOctahedron<ScalarT>;
    using HighSwRegularizationMethod = RegularizedTwoPLocalRulesBase::HighSwRegularizationMethod;

public:

    using Scalar = ScalarT;
    using Params = typename RegularizedTwoPLocalRulesBase::Params<Scalar>;

    static constexpr bool supportsMultipleGeometries()
    { return false; }

    /*!
     * \brief A regularized pore network capillary pressure-saturation
     *        curve.
     *
     * regularized part:
     *    - low saturation:  extend the \f$\mathrm{p_c(S_w)}\f$ curve with the slope at the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line (yes, there is a kink :-( ).
     *
     * For the non-regularized part:
     *
     * \copydetails PNMLocalRules::pc()
     */
    static Scalar pc(const Params& params, const Scalar sw)
    {
        const Scalar lowSw = params.lowSw;
        const Scalar highSw = params.highSw;

        // make sure that the capilary pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative is calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (sw < lowSw)
            return LocalRules::pc(params, lowSw) + mLow_(params) * (sw - lowSw);

        auto linearCurveForHighSw = [&]()
        {
            const Scalar slopeHighSw = -LocalRules::pc(params, highSw) / (1.0-highSw);
            return slopeHighSw*(sw - 1.0);
        };

        if (sw <= highSw)
            return LocalRules::pc(params, sw); // standard
        else if (sw <= 1.0) // regularized part below sw = 1.0
        {
            using std::pow;
            if (params.highSwRegularizationMethod == HighSwRegularizationMethod::powerLaw)
                return LocalRules::pc(params, highSw) * pow(((1.0-sw)/(1.0-highSw)), 1.0/3.0);

            else if (params.highSwRegularizationMethod == HighSwRegularizationMethod::linear)
                return linearCurveForHighSw();

            else if (params.highSwRegularizationMethod == HighSwRegularizationMethod::spline)
            {
                // use spline between threshold swe and 1.0
                const Scalar yTh = LocalRules::pc(params, highSw);
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

     /*! \brief The wetting-phase saturation of a pore body
     *
     * \copydetails PNMLocalRules::sw()
     */
    static Scalar sw(const Params& params, const Scalar pc)
    {
        // retrieve the low and the high threshold saturations for the
        // unregularized capillary pressure curve from the parameters
        const Scalar lowSw = params.lowSw;
        const Scalar highSw = params.highSw;
        const Scalar pcLowSw = LocalRules::pc(params, lowSw);
        const Scalar pcHighSw = LocalRules::pc(params, highSw);

        // low saturation, high pc:
        if (pc > pcLowSw)
            return lowSw + 1/mLow_(params)*(pc - pcLowSw);

        // high saturation, low pc:
        if (pc < pcHighSw)
        {
            if (params.highSwRegularizationMethod == HighSwRegularizationMethod::powerLaw)
            {
                auto result = pc/pcHighSw * pc/pcHighSw * pc/pcHighSw * (1.0-highSw);
                return 1.0 - result;
            }
            else if (params.highSwRegularizationMethod == HighSwRegularizationMethod::linear)
            {
                const Scalar slopeHighSw = -pcHighSw / (1.0-highSw);
                return pc/slopeHighSw + 1.0;
            }
            else if (params.highSwRegularizationMethod == HighSwRegularizationMethod::spline)
            {
                const Scalar yTh = LocalRules::pc(params, highSw);
                // invert spline between threshold swe and 1.0
                Spline<Scalar> sp(highSw, 1.0, // x0, x1
                                  yTh, 0, // y0, y1
                                  0.0, 0.0); // m0, m1

                return sp.intersectInterval(highSw, 1.0,
                                            0, 0, 0, pc);
            }
            else
                DUNE_THROW(Dune::NotImplemented, "Regularization not method not implemented");
        }

        return LocalRules::sw(params, pc);
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the wetting phase saturation.
     *
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    static Scalar dpc_dsw(const Params& params, const Scalar sw)
    {   // TODO!!!!
        assert(0 <= sw && sw <= 1);
        assert(params.shape == Pore::Shape::dodecahedron);
        using std::exp;
        const Scalar sigma = params.surfaceTension;
        const Scalar poreRadius = params.poreRadius;
        const Scalar e = exp(22.87*sw);
        return -(45.74*sigma*e) / (poreRadius*(e-1.0)*(e-1.0));
    }

    /*!
     * \brief DOCU
     *
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    static Scalar dsw_dpc(const Params& params, const Scalar sw)
    {
        return 0; // TODO
    }


private:

    /*!
     * \brief   The slope of the straight line used to regularize
     *          saturations below the minimum saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar mLow_(const Params& params)
    {
        return LocalRules::dpc_dsw(params, params.lowSw);
    }
};

}

#endif
