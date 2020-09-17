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
#ifndef REGULARIZED_PNM_2P_LOCAL_RULES_HH
#define REGULARIZED_PNM_2P_LOCAL_RULES_HH

#include "localrules.hh"

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
template<class Scalar, bool useZeroPc = true, class RegularizedLocalRulesForCube = TwoPLocalRulesCubeJoekarNiasar<Scalar>>
class RegularizedTwoPLocalRules : public TwoPLocalRulesBase
{

public:

    static constexpr bool supportsMultipleGeometries()
    { return true; }

    struct Params : public TwoPLocalRulesBase::Params<Scalar>
    {
        Scalar lowSw, highSw, slopeHighSw;
    };

    static Params makeParams(const Scalar poreRadius, const Scalar contactAngle, const Scalar surfaceTension, const Pore::Shape shape)
    {
        static const Scalar lowSw = getParam<Scalar>("Regularization.LowSw", 1e-2);
        static const Scalar highSw = getParam<Scalar>("Regularization.HighSw", 0.95);
        static const Scalar slopeHighSw = getParam<Scalar>("Regularization.SlopeHighSw", -1e9);
        return Params{poreRadius, contactAngle, surfaceTension, shape, lowSw, highSw, slopeHighSw};
    }

    /*!
     * \brief A regularized pore network capillary pressure-saturation
     *        curve.
     *
     * regularized part:
     *    - low saturation:  extend the \f$\mathrm{p_c(S_w)}\f$ curve with the slope at the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                       by a straight line (yes, there is a kink :-( ).
     *
     */
    static Scalar pc(const Params& params, const Scalar sw)
    {
        switch (params.shape)
        {
            case Pore::Shape::cube:
                return RegularizedLocalRulesForCube::pc(params, sw);
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
    }

     /*! \brief The wetting-phase saturation of a pore body
     *
     * \copydetails PNMLocalRules::sw()
     */
    static Scalar sw(const Params& params, const Scalar pc)
    {
        switch (params.shape)
        {
            case Pore::Shape::cube:
                return RegularizedLocalRulesForCube::sw(params, pc);
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
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
    {
        switch (params.shape)
        {
            case Pore::Shape::cube:
                return RegularizedLocalRulesForCube::dpc_dsw(params, sw);
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
    }
};

}

#endif
