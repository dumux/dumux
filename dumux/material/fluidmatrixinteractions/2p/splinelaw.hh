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
 * \brief Implementation of the capillary pressure and
 * relative permeability <-> saturation relations by spline interpolation.
 */
#ifndef DUMUX_SPLINE_LAW_HH
#define DUMUX_SPLINE_LAW_HH

#include "splinelawparams.hh"

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 *
 * \brief Implementation of the Brooks-Corey capillary pressure and
 * relative permeability <-> saturation relations by spline interpolation.
 *
 *\see SplineLawParams
 */
template <class ScalarT, class ParamsT = SplineLawParams<ScalarT> >
class SplineLaw
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    /*!
     * \brief The capillary pressure-saturation curve using spline interpolation.
     *
     * \param swe Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params Coefficients for the respective law.
     * \return Capillary pressure calculated by spline interpolation.
     */
    static Scalar pc(const Params &params, Scalar swe)
    {
        return params.pcSpline().eval(swe);
    }

    /*!
     * \brief The saturation-capillary pressure curve according to Brooks & Corey.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f$\mathrm{ \overline{S}_w = (\frac{p_C}{p_e})^{-\lambda}}\f$
     *
     * \param pc Capillary pressure \f$\mathrm{[p_C]}\f$  in \f$\mathrm{[Pa]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective wetting phase saturation calculated as inverse of SplineLaw constitutive relation.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar sw(const Params &params, Scalar pc)
    {
        DUNE_THROW(Dune::NotImplemented, "spline law sw-pc relation");
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar endPointPc(const Params &params)
    {
        return params.pcSpline().eval(1.0);
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation according to Brooks & Corey.
     *
     * This is equivalent to
     * \f$\mathrm{\frac{\partial p_C}{\partial \overline{S}_w} =
     * -\frac{p_e}{\lambda} \overline{S}_w^{-1/\lambda - 1}
     * }\f$
     *
     * \param swe Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Partial derivative of \f$\mathrm{[p_c]}\f$ w.r.t. effective saturation according to Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dpc_dswe(const Params &params, Scalar swe)
    {
        return params.pcSpline().evalDerivative(swe);
    }

    /*!
     * \brief The partial derivative of the effective
     *        saturation w.r.t. the capillary pressure according to Brooks & Corey.
     *
     * \param pc Capillary pressure \f$\mathrm{[p_c]}\f$ in \f$\mathrm{[Pa]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Partial derivative of effective saturation w.r.t. \f$\mathrm{[p_c]}\f$ according to Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dswe_dpc(const Params &params, Scalar pc)
    {
        DUNE_THROW(Dune::NotImplemented, "spline law dswe_dpc");
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by the Brooks-Corey
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Relative permeability of the wetting phase calculated as implied by Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar krw(const Params &params, Scalar swe)
    {
        return params.krwSpline().eval(swe);
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase with regard to the wetting saturation of the
     *        medium implied by the Brooks-Corey parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the relative permeability of the wetting phase w.r.t. effective wetting phase
     *                  saturation calculated as implied by Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dkrw_dswe(const Params &params, Scalar swe)
    {
        return params.krwSpline().evalDerivative(swe);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Relative permeability of the non-wetting phase calculated as implied by Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar krn(const Params &params, Scalar swe)
    {
        return params.krnSpline().eval(swe);
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the relative permeability of the non-wetting phase w.r.t. effective wetting phase
     *                  saturation calculated as implied by Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dkrn_dswe(const Params &params, Scalar swe)
    {
        return params.krnSpline().evalDerivative(swe);
    }

};

} // end namespace Dumux

#endif // SPLINE_LAW_HH
