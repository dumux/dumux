// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Implementation of the capillary pressure and
 * relative permeability <-> saturation relations according to Brooks and Corey.
 *
 */
#ifndef DUMUX_BROOKS_COREY_HH
#define DUMUX_BROOKS_COREY_HH

#include "brookscoreyparams.hh"

#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief Implementation of the Brooks-Corey capillary pressure <->
 *        saturation relation. This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vice versa.
 *
 * For general info: EffToAbsLaw
 *
 *\see BrooksCoreyParams
 */
template <class ScalarT, class ParamsT = BrooksCoreyParams<ScalarT> >
class BrooksCorey
{
public:
    typedef ParamsT     Params;
    typedef typename    Params::Scalar Scalar;

    /*!
     * \brief The capillary pressure-saturation curve according to Brooks & Corey.
     *
     * The Brooks-Corey empirical  capillary pressure <-> saturation
     * function is given by
     *
     *  \f$\mathrm{
        p_C = p_e\overline{S}_w^{-1/\lambda}
    *  }\f$
    *
     * \param swe Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
                        and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Capillary pressure calculated by Brooks & Corey constitutive relation.
     */
    static Scalar pc(const Params &params, Scalar swe)
    {
        assert(0 <= swe && swe <= 1);

        return params.pe()*pow(swe, -1.0/params.lambda());
    }

    /*!
     * \brief The saturation-capillary pressure curve according to Brooks & Corey.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f$\mathrm{
     \overline{S}_w = (\frac{p_C}{p_e})^{-\lambda}}\f$
     *
     * \param pc Capillary pressure \f$\mathrm{[p_C]}\f$  in \f$\mathrm{[Pa]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective wetting phase saturation calculated as inverse of BrooksCorey constitutive relation.
     */
    static Scalar sw(const Params &params, Scalar pc)
    {
        assert(pc >= 0);

        Scalar tmp = pow(pc/params.pe(), -params.lambda());
        return std::min(std::max(tmp, Scalar(0.0)), Scalar(1.0));
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation according to Brooks & Corey.
     *
     * This is equivalent to
     * \f$\mathrm{
     \frac{\partial p_C}{\partial \overline{S}_w} =
     -\frac{p_e}{\lambda} \overline{S}_w^{-1/\lambda - 1}
     }\f$
     *
     * \param swe Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Partial derivative of \f$\mathrm{[p_c]}\f$ w.r.t. effective saturation according to Brooks & Corey.
    */
    DUNE_DEPRECATED_MSG("dpc_dsw(const Params &params, Scalar swe) is deprecated. Use dpc_dswe(const Params &params, Scalar swe) instead.")
    static Scalar dpc_dsw(const Params &params, Scalar swe)
    {
        return dpc_dswe(params, swe);
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation according to Brooks & Corey.
     *
     * This is equivalent to
     * \f$\mathrm{
     \frac{\partial p_C}{\partial \overline{S}_w} =
     -\frac{p_e}{\lambda} \overline{S}_w^{-1/\lambda - 1}
     }\f$
     *
     * \param swe Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Partial derivative of \f$\mathrm{[p_c]}\f$ w.r.t. effective saturation according to Brooks & Corey.
    */
    static Scalar dpc_dswe(const Params &params, Scalar swe)
    {
        assert(0 <= swe && swe <= 1);

        return - params.pe()/params.lambda() * pow(swe, -1/params.lambda() - 1);
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
     */
    DUNE_DEPRECATED_MSG("dsw_dpc(const Params &params, Scalar pc) is deprecated. Use dswe_dpc(const Params &params, Scalar pc) instead.")
    static Scalar dsw_dpc(const Params &params, Scalar pc)
    {
        return dswe_dpc(params, pc);
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
     */
    static Scalar dswe_dpc(const Params &params, Scalar pc)
    {
        assert(pc >= 0);

        return -params.lambda()/params.pe() * pow(pc/params.pe(), - params.lambda() - 1);
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
     */
    static Scalar krw(const Params &params, Scalar swe)
    {
        assert(0 <= swe && swe <= 1);

        return pow(swe, 2.0/params.lambda() + 3);
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
     */
    DUNE_DEPRECATED_MSG("dkrw_dsw(const Params &params, Scalar swe) is deprecated. Use dkrw_dswe(const Params &params, Scalar swe) instead.")
    static Scalar dkrw_dsw(const Params &params, Scalar swe)
    {
        return dkrw_dswe(params, swe);
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
     */
    static Scalar dkrw_dswe(const Params &params, Scalar swe)
    {
        assert(0 <= swe && swe <= 1);

        return (2.0/params.lambda() + 3)*pow(swe, 2.0/params.lambda() + 2);
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
     */
    static Scalar krn(const Params &params, Scalar swe)
    {
        assert(0 <= swe && swe <= 1);

        Scalar exponent = 2.0/params.lambda() + 1;
        Scalar tmp = 1. - swe;
        return tmp*tmp*(1. - pow(swe, exponent));
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
     */
    DUNE_DEPRECATED_MSG("dkrn_dsw(const Params &params, Scalar swe) is deprecated. Use dkrn_dswe(const Params &params, Scalar swe) instead.")
    static Scalar dkrn_dsw(const Params &params, Scalar swe)
    {
        return dkrn_dswe(params, swe);
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
     */
    static Scalar dkrn_dswe(const Params &params, Scalar swe)
    {
        assert(0 <= swe && swe <= 1);

        return
            2.0*(swe - 1)*(
                1 +
                pow(swe, 2.0/params.lambda())*(
                    1.0/params.lambda() + 1.0/2 -
                    swe*(1.0/params.lambda() + 1.0/2)
                    )
                );
    }

};
}

#endif // BROOKS_COREY_HH
