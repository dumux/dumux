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
 * \brief Implementation of a regularized version of the Brooks-Corey
 *        capillary pressure / relative permeability  <-> saturation relation.
 */
#ifndef REGULARIZED_BROOKS_COREY_HH
#define REGULARIZED_BROOKS_COREY_HH

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

#include "brookscorey.hh"
#include "regularizedbrookscoreyparams.hh"

#include <dumux/common/spline.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the regularized  Brooks-Corey
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
 * \see BrooksCorey
 */
template <class ScalarT, class ParamsT = RegularizedBrooksCoreyParams<ScalarT> >
class RegularizedBrooksCorey
{
    using BrooksCorey = Dumux::BrooksCorey<ScalarT, ParamsT>;

public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    /*!
     * \brief A regularized Brooks-Corey capillary pressure-saturation
     *        curve.
     *
     * regularized part:
     *    - low saturation:  extend the \f$\mathrm{p_c(S_w)}\f$ curve with the slope at the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line (yes, there is a kink :-( ).
     *
     * For the non-regularized part:
     *
     * \copydetails BrooksCorey::pc()
     */
    static Scalar pc(const Params &params, Scalar swe)
    {
        const Scalar sThres = params.thresholdSw();

        // make sure that the capilary pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative is calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (swe <= sThres) {
            Scalar m = BrooksCorey::dpc_dswe(params, sThres);
            Scalar pcsweLow = BrooksCorey::pc(params, sThres);
            return pcsweLow + m*(swe - sThres);
        }
        else if (swe > 1.0) {
            Scalar m = BrooksCorey::dpc_dswe(params, 1.0);
            Scalar pcsweHigh = BrooksCorey::pc(params, 1.0);
            return pcsweHigh + m*(swe - 1.0);
        }

        // if the effective saturation is in an 'reasonable'
        // range, we use the real Brooks-Corey law...
        return BrooksCorey::pc(params, swe);
    }

    /*!
     * \brief   A regularized Brooks-Corey saturation-capillary pressure curve.
     *
     * regularized part:
     *    - low saturation:  extend the \f$\mathrm{p_c(S_w)}\f$ curve with the slope at the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line (yes, there is a kink :-( ).
     *
     *  The according quantities are obtained by exploiting theorem of intersecting lines.
     *
     * For the non-regularized part:
     *
     * \copydetails BrooksCorey::sw()
     */
    static Scalar sw(const Params &params, Scalar pc)
    {
        const Scalar sThres = params.thresholdSw();

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized version of
        // the Brooks-Corey law
        Scalar swe = BrooksCorey::sw(params, pc);

        // make sure that the capilary pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (swe <= sThres) {
            // invert the low saturation regularization of pc()
            Scalar m = BrooksCorey::dpc_dswe(params, sThres);
            Scalar pcsweLow = BrooksCorey::pc(params, sThres);
            return sThres + (pc - pcsweLow)/m;
        }
        else if (swe > 1.0) {
            Scalar m = BrooksCorey::dpc_dswe(params, 1.0);
            Scalar pcsweHigh = BrooksCorey::pc(params, 1.0);
            return 1.0 + (pc - pcsweHigh)/m;
        }

        return BrooksCorey::sw(params, pc);
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar endPointPc(const Params &params)
    { return params.pe(); }

    /*!
     * \brief A regularized version of the partial derivative
     *        of the \f$\mathrm{p_c(\overline{S}_w)}\f$ w.r.t. effective saturation
     *        according to Brooks & Corey.
     *
     * regularized part:
     *    - low saturation:  use the slope of the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line and use that slope (yes, there is a kink :-( ).
     *
     * For the non-regularized part:
     *
     * \copydetails BrooksCorey::dpc_dswe()
     */
    static Scalar dpc_dswe(const Params &params, Scalar swe)
    {
        const Scalar sThres = params.thresholdSw();

        // derivative of the regualarization
        if (swe <= sThres) {
            // calculate the slope of the straight line used in pc()
            Scalar m = BrooksCorey::dpc_dswe(params, sThres);
            return m;
        }
        else if (swe > 1.0) {
            // calculate the slope of the straight line used in pc()
            Scalar m = BrooksCorey::dpc_dswe(params, 1.0);
            return m;
        }

        return BrooksCorey::dpc_dswe(params, swe);
    }

    /*!
     * \brief A regularized version of the partial derivative
     *        of the \f$\mathrm{\overline{S}_w(p_c)}\f$ w.r.t. cap.pressure
     *        according to Brooks & Corey.
     *
     *  regularized part:
     *    - low saturation:  use the slope of the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line and use that slope (yes, there is a kink :-( ).
     *
     * For the non-regularized part:
     *
     * \copydetails BrooksCorey::dswe_dpc()
     */
    static Scalar dswe_dpc(const Params &params, Scalar pc)
    {
        const Scalar sThres = params.thresholdSw();

        //instead of return value = inf, return a very large number
        if (params.pe() == 0.0)
        {
            return 1e100;
        }

        // calculate the saturation which corresponds to the
        // saturation in the non-regularized version of the
        // Brooks-Corey law
        Scalar swe;
        if (pc < 0)
            swe = 1.5; // make sure we regularize below
        else
            swe = BrooksCorey::sw(params, pc);

        // derivative of the regularization
        if (swe <= sThres) {
            // calculate the slope of the straight line used in pc()
            Scalar m = BrooksCorey::dpc_dswe(params, sThres);
            return 1/m;
        }
        else if (swe > 1.0) {
            // calculate the slope of the straight line used in pc()
            Scalar m = BrooksCorey::dpc_dswe(params, 1.0);
            return 1/m;
        }
        return 1.0/BrooksCorey::dpc_dswe(params, swe);
    }

    /*!
     * \brief   Regularized version of the  relative permeability
     *          for the wetting phase of
     *          the medium implied by the Brooks-Corey
     *          parameterization.
     *
     *  regularized part:
     *    - below \f$\mathrm{\overline{S}_w =0}\f$:                  set relative permeability to zero
     *    - above \f$\mathrm{\overline{S}_w =1}\f$:                  set relative permeability to one
     *    - between \f$\mathrm{ 0.95 \leq \overline{S}_w \leq 1}\f$:  use a spline as interpolation
     *
     *  For not-regularized part:
        \copydetails BrooksCorey::krw()
     */
    static Scalar krw(const Params &params, Scalar swe)
    {
        if (swe <= 0.0)
            return 0.0;
        else if (swe >= 1.0)
            return 1.0;

        return BrooksCorey::krw(params, swe);
    }

    /*!
     * \brief A regularized version of the derivative of the relative
     *        permeability for the wetting phase in regard to the wetting
     *        saturation of the medium implied by the Brooks-Corey parameterization.
     *
     * \copydetails BrooksCorey::dkrw_dswe()
     */
    static Scalar dkrw_dswe(const Params &params, Scalar swe)
    {
        // derivative of the regularization
        // the slope is zero below sw=0.0 and above sw=1.0
        if (swe <= 0.0)
            return 0.0;
        else if (swe >= 1.0)
            return 0.0;

        return BrooksCorey::dkrw_dswe(params, swe);
    }

    /*!
     * \brief   Regularized version of the  relative permeability
     *          for the nonwetting phase of
     *          the medium implied by the Brooks-Corey
     *          parameterization.
     *
     * regularized part:
     *    - below \f$\mathrm{\overline{S}_w =0}\f$:                  set relative permeability to zero
     *    - above \f$\mathrm{\overline{S}_w =1}\f$:                  set relative permeability to one
     *    - for \f$\mathrm{0 \leq \overline{S}_w \leq 0.05}\f$:     use a spline as interpolation
     *
     * \copydetails BrooksCorey::krn()
     *
     */
    static Scalar krn(const Params &params, Scalar swe)
    {
        if (swe >= 1.0)
            return 0.0;
        else if (swe <= 0.0)
            return 1.0;

        return BrooksCorey::krn(params, swe);
    }

    /*!
     * \brief A regularized version of the derivative of the relative permeability
     *        for the nonwetting phase in regard to the wetting saturation of
     *        the medium as implied by the Brooks-Corey parameterization.
     *
     * \copydetails BrooksCorey::dkrn_dswe()
     */
    static Scalar dkrn_dswe(const Params &params, Scalar swe)
    {
        // derivative of the regularization
        // the slope is zero below sw=0.0 and above sw=1.0
        if (swe <= 0)
            return 0.0;
        else if (swe >= 1)
            return 0.0;

        return BrooksCorey::dkrn_dswe(params, swe);
    }
};
} // end namespace Dumux

#endif
