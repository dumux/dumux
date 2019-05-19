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
 * \brief   Implementation of the capillary pressure and
 *          relative permeability <-> saturation relations according to van Genuchten.
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_VAN_GENUCHTEN_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_VAN_GENUCHTEN_HH

#include <cmath>
#include <algorithm>

#include <dumux/common/parameters.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the van Genuchten capillary pressure <->
 *        saturation relation, and relative permeability.
 *
 * \note Capillary pressure model from van Genuchten (1980),
 *       relative permeability model from Mualem (1976)
 */
class VanGenuchten
{

public:
    /*!
     * \brief The parameter type
     * \tparam Scalar The scalar type
     * \note The van Genuchten laws are parameterized with three parameters: \f$\mathrm{n, m, \alpha}\f$.
     *       For the respective formulas check out the description of the free function.
     */
    template<class Scalar>
    struct Params
    {
        Scalar n, m, alpha, l;
    };

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    template<class Scalar = double>
    static Params<Scalar> makeParams(const std::string& paramGroup)
    {
        const auto vgn = getParamFromGroup<Scalar>(paramGroup, "Vgn");
        const auto vgm = 1.0 - 1.0/vgn;
        const auto vgAlpha = getParamFromGroup<Scalar>(paramGroup, "VgAlpha");
        const auto vgl = getParamFromGroup<Scalar>(paramGroup, "Vgl", 0.5);
        return {vgn, vgm, vgAlpha, vgl};
    }

    /*!
     * \brief The capillary pressure-saturation curve according to van Genuchten.
     *
     * Van Genuchten's empirical capillary pressure <-> saturation
     * function is given by
     * \f$\mathrm{
     p_C = (\overline{S}_w^{-1/m} - 1)^{1/n}/\alpha
     }\f$
     * \param swe Effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     * \note Instead of undefined behaviour if swe is not in the valid range, we return a valid number,
     *       by clamping the input. Note that pc(swe = 0.0) = inf, have a look at RegularizedVanGenuchten if this is a problem.
     */
    template<class Scalar>
    static Scalar pc(Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar pc = pow(pow(swe, -1.0/params.m) - 1, 1.0/params.n)/params.alpha;
        return pc;
    }

    /*!
     * \brief The saturation-capillary pressure curve according to van Genuchten.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f$\mathrm{
     \overline{S}_w = {p_C}^{-1} = ((\alpha p_C)^n + 1)^{-m}
     }\f$
     *
     * \param pc Capillary pressure \f$\mathrm{p_C}\f$ in \f$\mathrm{[Pa]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     * \return The effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       i.e. sw(pc < 0.0) = 0.0, by clamping the input to the physical bounds.
     */
    template<class Scalar>
    static Scalar sw(Scalar pc, const Params<Scalar>& params)
    {
        using std::pow;
        using std::max;

        pc = max(pc, 0.0); // the equation below is undefined for negative pcs

        const Scalar sw = pow(pow(params.alpha*pc, params.n) + 1, -params.m);
        return sw;
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar endPointPc(const Params<Scalar>& params)
    { return 0.0; }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation according to van Genuchten.
     *
     * This is equivalent to
     * \f$\mathrm{
     \frac{\partial p_C}{\partial \overline{S}_w} =
     -\frac{1}{\alpha} (\overline{S}_w^{-1/m} - 1)^{1/n - }
     \overline{S}_w^{-1/m} / \overline{S}_w / m
     }\f$
     *
     * \param swe Effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     * \note Instead of undefined behaviour if swe is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar dpc_dswe(Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar powSwe = pow(swe, -1/params.m);
        return - 1.0/params.alpha * pow(powSwe - 1, 1.0/params.n - 1)/params.n
                                  * powSwe/swe/params.m;
    }

    /*!
     * \brief The partial derivative of the effective
     *        saturation to the capillary pressure according to van Genuchten.
     *
     * \param pc Capillary pressure \f$\mathrm{p_C}\f$ in \f$\mathrm{[Pa]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar dswe_dpc(Scalar pc, const Params<Scalar>& params)
    {
        using std::pow;
        using std::max;

        pc = max(pc, 0.0); // the equation below is undefined for negative pcs

        const Scalar powAlphaPc = pow(params.alpha*pc, params.n);
        return -pow(powAlphaPc + 1, -params.m-1)*params.m*powAlphaPc/pc*params.n;
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by van Genuchten / Mualem
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar krw(Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar r = 1.0 - pow(1.0 - pow(swe, 1.0/params.m), params.m);
        return pow(swe, params.l)*r*r;
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase in regard to the wetting saturation of the
     *        medium implied by the van Genuchten parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar dkrw_dswe(Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar x = 1.0 - pow(swe, 1.0/params.m);
        const Scalar xToM = pow(x, params.m);
        return (1.0 - xToM)/pow(swe, params.l-1) * ( (1.0 - xToM)*params.l + 2*xToM*(1.0-x)/x );
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar krn(Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        return pow(1 - swe, params.l) * pow(1 - pow(swe, 1.0/params.m), 2*params.m);
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the van Genuchten
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar dkrn_dswe(Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const auto sne = 1.0 - swe;
        const auto x = 1.0 - pow(swe, 1.0/params.m);
        return -pow(sne, params.l-1.0) * pow(x, 2*params.m) - 1.0) * ( params.l*x + 2.0*sne/swe*(1.0 - x) );
    }
};

} // end namespace Dumux::FluidMatrix

#endif
