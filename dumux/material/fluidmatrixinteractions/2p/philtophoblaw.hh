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
 * \ingroup Fluidmatrixinteractions
 * \brief This material law takes a material law defined for effective
 *        saturations and converts it to a material law defined on
 *        absolute saturations. It is valid for hydrophobic materials and is
 *        called with the non-wetting phase saturation and then calls the
 *        material law defined for the wetting phase saturation.
 */
#ifndef DUMUX_PHIL_TO_PHOB_LAW_HH
#define DUMUX_PHIL_TO_PHOB_LAW_HH

#include <dune/common/exceptions.hh>

#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslawparams.hh>

namespace Dumux
{
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief This material law takes a material law defined for effective
 *        saturations and converts it to a material law defined on absolute
 *        saturations.
 *
 *        Additionally, this class converts the material law defined for the wetting phase so that it
 *        can be called with the non-wetting phase saturation. This is needed for hydrophobic materials.
 *
 *        The idea: "material laws" (like VanGenuchten or BrooksCorey) are defined for effective saturations.
 *        The numeric calculations however are performed with absolute saturations. The EffToAbsLaw class gets
 *        the "material laws" actually used as well as the corresponding parameter container as template arguments.
 *
 *        The desired function (pc, Sw... ) of the actually used "material laws" gets the non-wetting
 *        saturation from the model. Here, the wetting saturation is calculated from it: \f$\mathrm{S_w = 1 - S_n}\f$.
 *        Then the effective wetting saturations are calculated and handed to the material law.
 *
 *        The following definition shows the pc-Sw-relationship for the case in which the x-axis
 *        is \f$\mathrm{S_w (=S_{water})}\f$ and the y-axis is \f$\mathrm{p_c = p_{non_wetting}-p_{wetting}}\f$.
 *        \image html pc_Sw_1.png
 *
 *        But DumuX actually uses the following curve, because \f$\mathrm{p_c = p_{non-water}-p_{water}}\f$ for the y-axis is used.
 *        \image html pc_Sw_2.png
 *
 *        This approach makes sure that in the "material laws" only effective wetting saturations are considered,
 *        which makes sense, as these laws only deal with effective saturations. This also allows for changing
 *        the calculation of the effective saturations easily, as this is subject of discussion / may be problem specific.
 *
 *        Additionally, handing over effective saturations to the "material laws" in stead of them calculating effective
 *        saturations prevents accidentally "converting twice".
 *
 *        Additionally, the values calculated by the actual material laws are then modified: \f$\mathrm{pc_{phob} = - pc_{phil} (1-S_n)}\f$.
 *
 *        This boils down to:
 *        - the actual material laws (linear, VanGenuchten...) do not need to deal with any kind of conversion
 *        - the definition of the material law in the spatial parameters is not really intuitive, but using it is:
 *          Hand in values, get back values, do not deal with conversion.
 *
 *        CAREFUL: here w and n still stand for wetting and non-wetting. In the hydrophobic model w stands for
 *        water (non-wetting) and n for non-water (wetting). Hence, Swr is the wetting phase (gas) residual
 *        saturation and needs to be set accordingly in the spatial parameters.
 */

template <class EffLawT, class AbsParamsT = EffToAbsLawParams<typename EffLawT::Params> >
class PhilToPhobLaw : EffToAbsLaw<EffLawT>
{
    using EffLaw = EffLawT;

public:
    using Params = AbsParamsT;
    using Scalar = typename EffLaw::Scalar;
    using ParentType = EffToAbsLaw<EffLawT>;
    /*!
     * \brief The capillary pressure-saturation curve.
     *
     *
     * \param sn Absolute saturation of the non-wetting phase \f$\mathrm{\overline{S}_n}\f$. It is converted to the
     *                  effective saturation of the wetting phase and then handed over to the material law actually
     *                  used for calculation.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Capillary pressure \f$\mathrm{p_c}\f$ in \f$\mathrm{[Pa]}\f$  calculated by specific constitutive relation
     *                  (EffLaw e.g. Brooks & Corey, van Genuchten, linear...)*-1 to account for hydrophobic material.
     *
     */
    static Scalar pc(const Params &params, Scalar sn)
    {
        return -1.0 * EffLaw::pc(params, 1.0 - ParentType::swToSwe(params, sn));
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * \param pc Capillary pressure \f$\mathrm{p_c}\f$ in \f$\mathrm{[Pa]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     *\return Absolute wetting phase saturation calculated as inverse of
     *                  (EffLaw e.g. Brooks & Corey, van Genuchten, linear...) constitutive relation.
     *
     * \return The absolute saturation of the wetting phase \f$\mathrm{S_w}\f$
     */
    static Scalar sw(const Params &params, Scalar pc)
    {
       DUNE_THROW(Dune::NotImplemented, "PhilToPhobLaw::sw(params, pc)");
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure w.r.t the absolute saturation.
     *
     *        In this case the chain rule needs to be applied:
     \f$\mathrm{
             p_c = p_c( \overline{S}_w (S_w))
             \rightarrow p_c ^\prime = \frac{\partial  p_c}{\partial \overline{S}_w} \frac{\partial \overline{S}_w}{\partial S_w}
     }\f$
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and
     *                  then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return          Partial derivative of \f$\mathrm{p_c}\f$ w.r.t. effective saturation according to EffLaw
     *                  e.g. Brooks & Corey, van Genuchten, linear... .
     */
    static Scalar dpc_dsw(const Params &params, Scalar sw)
    {
        DUNE_THROW(Dune::NotImplemented, "PhilToPhobLaw::dpc_dsw(params, sw)");
    }

    /*!
     * \brief Returns the partial derivative of the absolute
     *        saturation w.r.t. the capillary pressure.
     *
     * In this case the chain rule needs to be applied:
     \f$\mathrm{
            S_w = S_w(\overline{S}_w (p_c) )
            \rightarrow S_w^\prime = \frac{\partial S_w}{\partial \overline{S}_w} \frac{\partial \overline{S}_w}{\partial p_c}
    }\f$
     *
     *
     * \param pc Capillary pressure \f$\mathrm{p_c}\f$ in \f$\mathrm{[Pa]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Partial derivative of effective saturation w.r.t. \f$\mathrm{p_C}\f$ according to EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     */
    static Scalar dsw_dpc(const Params &params, Scalar pc)
    {
        DUNE_THROW(Dune::NotImplemented, "PhilToPhobLaw::dsw_dpc(params, pc)");
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param sn Absolute saturation of the non-wetting phase \f$\mathrm{\overline{S}_w}\f$. It is converted to effective saturation of the wetting phase
     *                  and then handed over to the material law actually used for calculation.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Relative permeability of the wetting phase calculated as implied by EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     *
     */
    static Scalar krw(const Params &params, Scalar sn)
    {
        if(sn >= params.swr())
            return sn*sn;//EffLaw::krn(params, SwToSwe(params, 1.0 - sn));
        else
            return 0.0;
    };

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     * \param sn Absolute saturation of the non-wetting phase \f$\mathrm{\overline{S}_w}\f$. It is converted to effective saturation of the wetting phase
     *                  and then handed over to the material law actually used for calculation.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Relative permeability of the non-wetting phase calculated as implied by EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     */
    static Scalar krn(const Params &params, Scalar sn)
    {
        if((1.0-sn) >= params.snr())
            return (1.0-sn)*(1.0-sn)*(1.0-sn);//EffLaw::krw(params, SwToSwe(params, 1.0 - sn));
        else
            return 0.0;
    }
};
}

#endif
