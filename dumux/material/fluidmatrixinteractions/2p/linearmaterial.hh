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
 * \brief   Linear capillary pressure and
 *          relative permeability <-> saturation relations
 */
#ifndef LINEAR_MATERIAL_HH
#define LINEAR_MATERIAL_HH

#include "linearmaterialparams.hh"

#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Linear capillary pressure and
 * relative permeability <-> saturation relations
 *
 * The entry pressure is reached at \f$\mathrm{ \overline{S}_w = 1}\f$, the maximum
 * capillary pressure is observed at \f$\mathrm{ \overline{S}_w = 0}\f$.
 *
 * For general info: EffToAbsLaw
 *
 * \see LinearMaterialParams
 */
template <class ScalarT, class ParamsT = LinearMaterialParams<ScalarT> >
class LinearMaterial
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    /*!
     * \brief The linear capillary pressure-saturation curve.
     *
     * This material law is linear:
     * \f$\mathrm{
     p_C = (1 - \overline{S}_w ) (p_{C,max} - p_{C,entry}) + p_{C,entry}
     }\f$
     *
     * \param swe Effective saturation of the wetting phase \f$\overline{S}_w\f$ conversion from absolute saturation happened in EffToAbsLaw.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Capillary pressure calculated by linear constitutive relation.
     */
    static Scalar pc(const Params &params, Scalar swe)
    {
        return (1 - swe)*(params.maxPc() - params.entryPc()) + params.entryPc();
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f$\mathrm{
         S_w = 1 - \frac{p_C - p_{C,entry}}{p_{C,max} - p_{C,entry}}
     }\f$
     *
     * \param pc Capillary pressure \f$\mathrm{[p_C]}\f$ in \f$\mathrm{[Pa]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective wetting phase saturation calculated as inverse of the linear constitutive relation.
     */
    static Scalar sw(const Params &params, Scalar pc)
    {
        return 1 - (pc - params.entryPc())/(params.maxPc() - params.entryPc());
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar endPointPc(const Params &params)
    { return params.entryPc(); }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure w.r.t. the effective saturation.
     *
     * This is equivalent to
     * \f$\mathrm{
     \frac{\partial p_C}{\partial \overline{S}_w} =
     - (p_{C,max} - p_{C,min})
     }\f$
     * \param swe  Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$ conversion from absolute saturation happened in EffToAbsLaw.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Partial derivative of \f$\mathrm{[p_c]}\f$ w.r.t. effective saturation according to linear material relation.
     */
    static Scalar dpc_dswe(const Params &params, Scalar swe)
    {
        return - (params.maxPc() - params.entryPc());
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     *
     * \param pc Capillary pressure \f$\mathrm{[p_C]}\f$  in \f$\mathrm{[Pa]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Partial derivative of effective saturation w.r.t. \f$\mathrm{[p_c]}\f$ according to linear relation.
     */
    static Scalar dsw_dpc(const Params &params, Scalar pc)
    {
        return - 1/(params.maxPc() - params.entryPc());
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \param swe Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$ conversion from absolute saturation happened in EffToAbsLaw.
     * \return Relative permeability of the wetting phase calculated as linear relation.
     */
    static Scalar krw(const Params &params, Scalar swe)
    {
        using std::max;
        using std::min;
        return max(min(swe,1.0),0.0);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \param swe Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$ conversion from absolute saturation happened in EffToAbsLaw.
     * \return Relative permeability of the non-wetting phase calculated as linear relation.
     */
    static Scalar krn(const Params &params, Scalar swe)
    {
        Scalar sne = 1 - swe;
        using std::max;
        using std::min;
        return max(min(sne,1.0),0.0);
    }
};
}

#endif
