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
 * \brief This material law takes a material law defined for effective
 *        saturations and converts it to a material law defined on
 *        absolute saturations.
 */
#ifndef DUMUX_EFF_TO_ABS_LAW_HH
#define DUMUX_EFF_TO_ABS_LAW_HH

#include <dune/common/exceptions.hh>
#include "efftoabslawparams.hh"

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief This material law takes a material law defined for effective
 *        saturations and converts it to a material law defined on absolute
 *        saturations.
 *
 *        The idea: "material laws" (like VanGenuchten or BrooksCorey) are defined for effective saturations.
 *        The numeric calculations however are performed with absolute saturations. The EffToAbsLaw class gets
 *        the "material laws" actually used as well as the corresponding parameter container as template arguments.
 *
 *        Subsequently, the desired function (pc, sw... ) of the actually used "material laws" are called but with the
 *        saturations already converted from absolute to effective.
 *
 *        This approach makes sure that in the "material laws" only effective saturations are considered, which makes sense,
 *        as these laws only deal with effective saturations. This also allows for changing the calculation of the effective
 *        saturations easily, as this is subject of discussion / may be problem specific.
 *
 *        Additionally, handing over effective saturations to the "material laws" in stead of them calculating effective
 *        saturations prevents accidently "converting twice".
 *
 *        This boils down to:
 *        - the actual material laws (linear, VanGenuchten...) do not need to deal with any kind of conversion
 *        - the definition of the material law in the spatial parameters is not really intuitive, but using it is:
 *          Hand in values, get back values, do not deal with conversion.
 */
template <class EffLawT, class AbsParamsT = EffToAbsLawParams<typename EffLawT::Params> >
class EffToAbsLaw
{
    using EffLaw = EffLawT;

public:
    using Params = AbsParamsT;
    using Scalar = typename EffLaw::Scalar;

    /*!
     * \brief The capillary pressure-saturation curve.
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$. It is converted to effective saturation
     *                  and then handed over to the material law actually used for calculation.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return          Capillary pressure calculated by specific constitutive relation
     *                  (EffLaw e.g. Brooks & Corey, van Genuchten, linear...)
     */
    static Scalar pc(const Params &params, const Scalar sw)
    {
        return EffLaw::pc(params, swToSwe(params, sw));
    }

    /*!
     * \brief The capillary pressure-saturation curve for the gas and wetting phase
     * \param params Array of parameters
     * \param sw wetting phase saturation or sum of wetting phase saturations
     */
    static Scalar pcgw(const Params &params, const Scalar sw)
    {
         return EffLaw::pcgw(params, swToSwe(params, sw));
    }

    /*!
     * \brief The capillary pressure-saturation curve the nonwetting and wetting phase
     * \param params Array of parameters
     * \param sw wetting phase saturation or sum of wetting phase saturations
     */
    static Scalar pcnw(const Params &params, const Scalar sw)
    {
        return EffLaw::pcnw(params, swToSwe(params, sw));
    }

    /*!
     * \brief The capillary pressure-saturation curve for the gas and nonwetting phase
     * \param params Array of parameters
     * \param st sum of wetting (liquid) phase saturations
     */
    static Scalar pcgn(const Params &params, const Scalar st)
    {
        return EffLaw::pcgn(params, stToSte(params, st));
    }

    /*!
     * \brief This function ensures a continuous transition from 2 to 3 phases and vice versa
     * \param params Array of parameters
     * \param sn Nonwetting liquid saturation
     */
    static Scalar pcAlpha(const Params &params, const Scalar sn)
    {
        return EffLaw::pcAlpha(params, sn);
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     * \param pc Capillary pressure \f$\mathrm{[p_c]}\f$ in \f$\mathrm{[Pa]}\f$.
     * \param params  A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     *\return           Absolute wetting phase saturation calculated as inverse of
     *                  (EffLaw e.g. Brooks & Corey, van Genuchten, linear...) constitutive relation.
     *
     * \return The absolute saturation of the wetting phase \f$\mathrm{[S_w]}\f$
     */
    static Scalar sw(const Params &params, const Scalar pc)
    {
        return EffLaw::sw(params, pc);
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure w.r.t the absolute saturation.
     *        In this case the chain rule needs to be applied:
     \f$\mathrm{
             p_c = p_c( \overline{S}_w (S_w))
             \rightarrow p_c ^\prime = \frac{\partial  p_c}{\partial \overline S_w} \frac{\partial \overline{S}_w}{\partial S_w}
     }\f$
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
                        and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return          Partial derivative of \f$\mathrm{[p_c]}\f$ w.r.t. effective saturation according to
                        EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     */
    static Scalar dpc_dsw(const Params &params, const Scalar sw)
    {
        return EffLaw::dpc_dswe(params, pc);
    }

    /*!
     * \brief Returns the partial derivative of the absolute
     *        saturation w.r.t. the capillary pressure.
     * In this case the chain rule needs to be applied:
     \f$\mathrm{
            S_w = S_w(\overline{S}_w (p_c) )
            \rightarrow S_w^\prime = \frac{\partial S_w}{\partial \overline{S}_w} \frac{\partial \overline{S}_w}{\partial p_c}
     }\f$
     *
     *
     * \param pc Capillary pressure \f$\mathrm{[p_c]}\f$ in \f$\mathrm{[Pa]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
                        and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Partial derivative of effective saturation w.r.t. \f$\mathrm{[p_c]}\f$ according to
                        EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     */
    static Scalar dsw_dpc(const Params &params, const Scalar pc)
    {
        return EffLaw::dsw_dpc(params, pc);
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$. It is converted to effective saturation
     *                  and then handed over to the material law actually used for calculation.
     * \param sn Absolute saturation of the nonwetting phase \f$\mathrm{[{S}_n]}\f$. It is converted to effective saturation
     *                  and then handed over to the material law actually used for calculation.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Relative permeability of the wetting phase calculated as implied by
     *                  EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     *
     */
    static Scalar krw(const Params &params, Scalar sw, const Scalar sn)
    {
        return EffLaw::krw(params, swToSwe(params, sw));
    }

    /*!
     * \brief The relative permeability for the nonwetting phase.
     *
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{[{S}_w]}\f$. It is converted to effective saturation
     *                  and then handed over to the material law actually used for calculation.
     * \param sn Absolute saturation of the nonwetting phase \f$\mathrm{[{S}_n]}\f$. It is converted to effective saturation
     *                  and then handed over to the material law actually used for calculation.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return          Relative permeability of the nonwetting phase calculated as implied by
     *                  EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     */
    static Scalar krn(const Params &params, const Scalar sw, const Scalar sn)
    {
        const Scalar st = sw+sn;
        return EffLaw::krn(params, swToSwe(params, sw), sn, stToSte(params, st));
    }

    /*!
     * \brief The relative permeability for the gas phase.
     *
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{[{S}_w]}\f$. It is converted to effective saturation
     *                  and then handed over to the material law actually used for calculation.
     * \param sn Absolute saturation of the nonwetting phase \f$\mathrm{[{S}_n]}\f$. It is converted to effective saturation
     *                  and then handed over to the material law actually used for calculation.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return          Relative permeability of the nonwetting phase calculated as implied by
     *                  EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     */
    static Scalar krg(const Params &params, const Scalar sw, const Scalar sn)
    {
        const Scalar st = sw+sn;
        return EffLaw::krg(params, stToSte(params, st));
    }

    /*!
     * \brief The relative permeability for a phase.
     * \param sw Wetting liquid saturation
     * \param sg Gas saturation
     * \param sn Nonwetting liquid saturation
     * \param params Array of parameters.
     * \param phaseIdx indicator, The saturation of all phases.
     */
    static Scalar kr(const Params &params, const int phaseIdx, const Scalar sw, const Scalar sn, const Scalar sg)
    {
      const Scalar st = sw+sn;
      return EffLaw::kr(params, phaseIdx, swToSwe(params, sw), sn, stToSte(params, st));
    }

    /*!
     * \brief the basis for calculating adsorbed NAPL in storage term
     * \param params Array of parameters
     */
    static Scalar bulkDensTimesAdsorpCoeff(const Params &params)
    {
        return EffLaw::bulkDensTimesAdsorpCoeff(params);
    }

    /*!
     * \brief Convert an absolute wetting saturation to an effective one.
     *
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{[{S}_w]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective saturation of the wetting phase.
     */
    static Scalar swToSwe(const Params &params, const Scalar sw)
    {
       return (sw-params.swr())/(1.-params.swr());
    }

    /*!
     * \brief Convert an absolute nonwetting saturation to an effective one.
     *
     * \param sn Absolute saturation of the nonwetting phase \f$\mathrm{[{S}_n]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective saturation of the nonwetting phase.
     */
    static Scalar snToSne(const Params &params, const Scalar sn)
    {
        return sn; // sne equals sn
    }

    /*!
     * \brief Convert an absolute total liquid saturation to an effective one.
     *
     * \param st Absolute saturation of the total liquid phase (sw+sn) \f$\mathrm{[{S}_n]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective saturation of the nonwetting phase.
     */
    static Scalar stToSte(const Params &params, const Scalar st)
    {
        return (st-params.swr()) / (1-params.swr());
    }

    /*!
     * \brief Convert an absolute gas saturation to an effective one.
     *
     * \param sg Absolute saturation of the gas phase \f$\mathrm{[{S}_n]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective saturation of the nonwetting phase.
     */
    static Scalar sgToSge(const Params &params, Scalar sg)
    {
        DUNE_THROW(Dune::NotImplemented, "sgTosge for three phases not implemented!");
    }

//private:
    /*!
     * \brief Convert an effective wetting saturation to an absolute one.
     *
     * \param swe Effective saturation of the nonwetting phase \f$\mathrm{[\overline{S}_n]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Absolute saturation of the nonwetting phase.
     */
    static Scalar sweToSw_(const Params &params, Scalar swe)
    {
         DUNE_THROW(Dune::NotImplemented, "sweTosw for three phases not implemented!");
    }

    static Scalar sneToSn_(const Params &params, Scalar sne)
    {
         DUNE_THROW(Dune::NotImplemented, "sneTosn for three phases not implemented!");
    }

    static Scalar sgeToSg_(const Params &params, Scalar sge)
    {
         DUNE_THROW(Dune::NotImplemented, "sgeTosg for three phases not implemented!");
    }
    /*!
     * \brief Derivative of the effective saturation w.r.t. the absolute saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the effective saturation w.r.t. the absolute saturation.
     */
    static Scalar dswe_dsw_(const Params &params)
    {
       DUNE_THROW(Dune::NotImplemented, "dswe_dsw for three phases not implemented!");
    }

    /*!
     * \brief Derivative of the absolute saturation w.r.t. the effective saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the absolute saturation w.r.t. the effective saturation.
     */
    static Scalar dsw_dswe_(const Params &params)
    {
       DUNE_THROW(Dune::NotImplemented, "dsw_dswe for three phases not implemented!");
    }
};
} // end namespace Dumux

#endif
