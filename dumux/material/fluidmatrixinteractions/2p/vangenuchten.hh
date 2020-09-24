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

// remove from here after release 3.3 /////////////
#include "vangenuchtenparams.hh"

#include <algorithm>
#include <cmath>
#include <cassert>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the van Genuchten capillary pressure <->
 *        saturation relation. This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vice versa.
 *
 * For general info: EffToAbsLaw
 *
 * \note Capillary pressure model from van Genuchten (1980),
 *       relative permeability model from Mualem (1976)
 *
 * \see VanGenuchtenParams
 */
template <class ScalarT, class ParamsT = VanGenuchtenParams<ScalarT> >
class [[deprecated("Use new material laws and FluidMatrix::VanGenuchten instead!")]] VanGenuchten
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

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
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \note Instead of undefined behaviour if swe is not in the valid range, we return a valid number,
     *       by clamping the input. Note that for pc(swe = 0.0) = inf, have a look at RegularizedVanGenuchten if this is a problem.
     */
    static Scalar pc(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar pc = pow(pow(swe, -1.0/params.vgm()) - 1, 1.0/params.vgn())/params.vgAlpha();
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
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          The effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       i.e. sw(pc < 0.0) = 0.0, by clamping the input to the physical bounds.
     */
    static Scalar sw(const Params &params, Scalar pc)
    {
        using std::pow;
        using std::max;

        pc = max(pc, 0.0); // the equation below is undefined for negative pcs

        const Scalar sw = pow(pow(params.vgAlpha()*pc, params.vgn()) + 1, -params.vgm());
        return sw;
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first,
     *                  the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar endPointPc(const Params &params)
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
     *                  Therefore, in the (problem specific) spatialParameters
     *                  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if swe is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dpc_dswe(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar powSwe = pow(swe, -1/params.vgm());
        return - 1.0/params.vgAlpha() * pow(powSwe - 1, 1.0/params.vgn() - 1)/params.vgn()
                                      * powSwe/swe/params.vgm();
    }

    /*!
     * \brief The partial derivative of the effective
     *        saturation to the capillary pressure according to van Genuchten.
     *
     * \param pc Capillary pressure \f$\mathrm{p_C}\f$ in \f$\mathrm{[Pa]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters
     *                  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dswe_dpc(const Params &params, Scalar pc)
    {
        using std::pow;
        using std::max;

        pc = max(pc, 0.0); // the equation below is undefined for negative pcs

        const Scalar powAlphaPc = pow(params.vgAlpha()*pc, params.vgn());
        return -pow(powAlphaPc + 1, -params.vgm()-1)*params.vgm()*powAlphaPc/pc*params.vgn();
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by van Genuchten / Mualem
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters
     *                  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar krw(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar r = 1.0 - pow(1.0 - pow(swe, 1.0/params.vgm()), params.vgm());
        return pow(swe, params.vgl())*r*r;
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase in regard to the wetting saturation of the
     *        medium implied by the van Genuchten / Mualem parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters
     *                  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dkrw_dswe(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar x = 1.0 - pow(swe, 1.0/params.vgm());
        const Scalar xToM = pow(x, params.vgm());
        return (1.0 - xToM)*pow(swe, params.vgl()-1) * ( (1.0 - xToM)*params.vgl() + 2*xToM*(1.0-x)/x );
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters
     *                  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     * \note See e.g. Dury, Fischer, Schulin (1999) for application of Mualem model to non-wetting rel. perm.
     */
    static Scalar krn(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        return pow(1 - swe, params.vgl()) * pow(1 - pow(swe, 1.0/params.vgm()), 2*params.vgm());
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the van Genuchten
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters
     *                  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dkrn_dswe(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const auto sne = 1.0 - swe;
        const auto x = 1.0 - pow(swe, 1.0/params.vgm());
        return -pow(sne, params.vgl()-1.0) * pow(x, 2*params.vgm() - 1.0) * ( params.vgl()*x + 2.0*sne/swe*(1.0 - x) );
    }

};

} // end namespace Dumux
// remove until here after release 3.3 /////////////

#include <cmath>
#include <algorithm>

#include <dumux/common/parameters.hh>
#include <dumux/common/spline.hh>
#include <dumux/common/optionalscalar.hh>
#include <dumux/material/fluidmatrixinteractions/2p/materiallaw.hh>

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
     * \note The van Genuchten laws are parameterized with four parameters: \f$\mathrm{n, m, \alpha, l}\f$.
     *
     * - \f$\mathrm{\alpha}\f$ shape parameter \f$\mathrm{[1/Pa]}\f$
     * - \f$\mathrm{m}\f$ shape parameter \f$\mathrm{[-]}\f$
     * - \f$\mathrm{n}\f$ shape parameter \f$\mathrm{[-]}\f$
     * - \f$\mathrm{l}\f$ pore-connectivity parameter \f$\mathrm{[-]}\f$ of Mualem's relative permeability curve
     *
     * \note In the orignal Mualem (1976) paper the pore-connectivity parameter is called "n". It's referred to as "l" in
     *       several later publication of van Genuchten, e.g. van Genuchten (1991), Shaap & van Genuchten (2006).
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
        // l is usually chosen to be 0.5 (according to Mualem (1976), WRR)
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
    static Scalar swe(Scalar pc, const Params<Scalar>& params)
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
        return (1.0 - xToM)*pow(swe, params.l-1) * ( (1.0 - xToM)*params.l + 2*xToM*(1.0-x)/x );
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
        return -pow(sne, params.l-1.0) * pow(x, 2*params.m - 1.0) * ( params.l*x + 2.0*sne/swe*(1.0 - x) );
    }
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A regularization for the VanGenuchten material law
 * \note Regularization can be turned of by setting the threshold parameters
 *       out of range (runtime) or by replacing
 *       this class by NoRegularization (compile time).
 */
template <class Scalar>
class VanGenuchtenRegularization
{
public:
    //! Regularization parameters
    template<class S>
    struct Params
    {
        S pcLowSw, pcHighSw;
        S krwHighSw, krnLowSw;
    };

    //! Initialize the spline
    template<class MaterialLaw>
    void init(const MaterialLaw* m, const std::string& paramGroup)
    {
        pcLowSw_ = getParamFromGroup<Scalar>(paramGroup, "VgPcLowSwThreshold", 0.01);
        pcHighSw_ = getParamFromGroup<Scalar>(paramGroup, "VgPcHighSwThreshold", 0.99);
        krwHighSw_ = getParamFromGroup<Scalar>(paramGroup, "VgKrwHighSwThreshold", 0.1);
        krnLowSw_ = getParamFromGroup<Scalar>(paramGroup, "VgKrnLowSwThreshold", 0.9);

        initPcParameters_(m, pcLowSw_, pcHighSw_);
        initKrParameters_(m, krnLowSw_, krwHighSw_);
    }

    template<class MaterialLaw, class BaseParams, class EffToAbsParams>
    void init(const MaterialLaw* m, const BaseParams& bp, const EffToAbsParams& etap, const Params<Scalar>& p)
    {
        pcLowSw_ = p.pcLowSw;
        pcHighSw_ = p.pcHighSw;
        krwHighSw_ = p.krwHighSw;
        krnLowSw_ = p.krnLowSw;

        initPcParameters_(m, pcLowSw_, pcHighSw_);
        initKrParameters_(m, krnLowSw_, krwHighSw_);
    }

    /*!
     * \brief The regularized capillary pressure-saturation curve
     * regularized part:
     *  - low saturation:  extend the \f$\mathrm{p_c(S_w)}\f$ curve with the slope at the regularization point (i.e. no kink).
     *  - high saturation: connect the high regularization point with \f$\mathrm{\overline{S}_w =1}\f$
     *                     with a spline and continue linearly for \f$\mathrm{\overline{S}_w > 1}\f$
     */
    OptionalScalar<Scalar> pc(const Scalar swe) const
    {
        // make sure that the capillary pressure observes a derivative
        // != 0 for 'illegal' saturations. This is favourable for the
        // newton solver (if the derivative is calculated numerically)
        // in order to get the saturation moving to the right
        // direction if it temporarily is in an 'illegal' range.
        if (swe < pcLowSw_)
            return pcLowSwPcValue_ + pcDerivativeLowSw_*(swe - pcLowSw_);

        else if (swe > 1.0)
            return pcDerivativeHighSwEnd_*(swe - 1.0);

        else if (swe > pcHighSw_)
            return pcSpline_.eval(swe);

        else
            return {}; // no regularization
    }

    /*!
     * \brief The regularized partial derivative of the capillary pressure w.r.t. the saturation
     */
    OptionalScalar<Scalar> dpc_dswe(const Scalar swe) const
    {
        if (swe < pcLowSw_)
            return pcDerivativeLowSw_;

        else if (swe > 1.0)
            return pcDerivativeHighSwEnd_;

        else if (swe > pcHighSw_)
            return pcSpline_.evalDerivative(swe);

        else
            return {}; // no regularization
    }

    /*!
     * \brief The regularized saturation-capillary pressure curve
     */
    OptionalScalar<Scalar> swe(const Scalar pc) const
    {
        if (pc <= 0.0)
        {
            if (pcHighSw_ > 1.0 - std::numeric_limits<Scalar>::epsilon())
                return 1.0;
            else
                return pc/pcDerivativeHighSwEnd_ + 1.0;
        }

        // invert spline
        else if (pc <  pcHighSwPcValue_)
            return pcSpline_.intersectInterval(pcHighSw_, 1.0, 0.0, 0.0, 0.0, pc);

        else if (pc >= pcLowSwPcValue_)
            return (pc - pcLowSwPcValue_)/pcDerivativeLowSw_ + pcLowSw_;

        else
            return {}; // no regularization
    }

    /*!
     * \brief The regularized partial derivative of the saturation to the capillary pressure
     */
    OptionalScalar<Scalar> dswe_dpc(const Scalar pc) const
    {
        if (pc <= 0.0)
        {
            if (pcHighSw_ > 1.0 - std::numeric_limits<Scalar>::epsilon())
                return 0.0;
            else
                return 1.0/pcDerivativeHighSwEnd_;
        }

        // derivative of the inverse of the function is one over derivative of the function
        else if (pc <  pcHighSwPcValue_)
            return 1.0/pcSpline_.evalDerivative(pcSpline_.intersectInterval(pcHighSw_, 1.0, 0.0, 0.0, 0.0, pc));

        else if (pc >= pcLowSwPcValue_)
            return 1.0/pcDerivativeLowSw_;

        else
            return {}; // no regularization
    }

    /*!
     * \brief The regularized relative permeability for the wetting phase
     */
    OptionalScalar<Scalar> krw(const Scalar swe) const
    {
        if (swe < 0.0)
            return 0.0;
        else if (swe > 1.0 - std::numeric_limits<Scalar>::epsilon())
            return 1.0;
        else if (swe > krwHighSw_)
            return krwSpline_.eval(swe);
        else
            return {}; // no regularization
    }

    /*!
     * \brief The regularized derivative of the relative permeability for the wetting phase w.r.t. saturation
     */
    OptionalScalar<Scalar> dkrw_dswe(const Scalar swe) const
    {
        if (swe < 0.0)
            return 0.0;
        else if (swe > 1.0 - std::numeric_limits<Scalar>::epsilon())
            return 0.0;
        else if (swe > krwHighSw_)
            return krwSpline_.evalDerivative(swe);
        else
            return {}; // no regularization
    }

    /*!
     * \brief The regularized relative permeability for the non-wetting phase
     */
    OptionalScalar<Scalar> krn(const Scalar swe) const
    {
        if (swe < 0.0)
            return 1.0;
        else if (swe > 1.0 - std::numeric_limits<Scalar>::epsilon())
            return 0.0;
        else if (swe < krnLowSw_)
            return krnSpline_.eval(swe);
        else
            return {}; // no regularization
    }

    /*!
     * \brief The regularized derivative of the relative permeability for the non-wetting phase w.r.t. saturation
     */
    OptionalScalar<Scalar> dkrn_dswe(const Scalar swe) const
    {
        if (swe < 0.0)
            return 0.0;
        else if (swe > 1.0 - std::numeric_limits<Scalar>::epsilon())
            return 0.0;
        else if (swe < krnLowSw_)
            return krnSpline_.evalDerivative(swe);
        else
            return {}; // no regularization
    }

private:
    template<class MaterialLaw>
    void initPcParameters_(const MaterialLaw* m, const Scalar lowSw, const Scalar highSw)
    {
        pcDerivativeLowSw_ = m->template dpc_dswe<false>(lowSw);

        pcDerivativeHighSwThreshold_ = m->template dpc_dswe<false>(highSw);
        pcDerivativeHighSwEnd_ = 2.0*(0.0 - m->template pc<false>(highSw))/(1.0 - highSw);

        pcLowSwPcValue_ = m->template pc<false>(lowSw);
        pcHighSwPcValue_ = m->template pc<false>(highSw);

        pcSpline_ = Spline<Scalar>(highSw, 1.0, // x0, x1
                                   pcHighSwPcValue_, 0, // y0, y1
                                   pcDerivativeHighSwThreshold_, pcDerivativeHighSwEnd_); // m0, m1


    }

    template<class MaterialLaw>
    void initKrParameters_(const MaterialLaw* m, const Scalar lowSw, const Scalar highSw)
    {
        krwSpline_ = Spline<Scalar>(highSw, 1.0, // x0, x1
                                    m->template krw<false>(highSw), 1.0, // y0, y1
                                    m->template dkrw_dswe<false>(highSw), 0.0); // m0, m1

        krnSpline_ = Spline<Scalar>(0.0, lowSw, // x0, x1
                                    1.0, m->template krn<false>(lowSw), // y0, y1
                                    0.0, m->template dkrn_dswe<false>(lowSw)); // m0, m1

    }

    Scalar pcLowSw_, pcHighSw_;
    Scalar pcLowSwPcValue_, pcHighSwPcValue_;
    Scalar krwHighSw_, krnLowSw_;
    Scalar pcDerivativeLowSw_;
    Scalar pcDerivativeHighSwThreshold_, pcDerivativeHighSwEnd_;

    Spline<Scalar> pcSpline_;
    Spline<Scalar> krwSpline_;
    Spline<Scalar> krnSpline_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A default configuration for using the VanGenuchten material law
 */
template<typename Scalar = double>
using VanGenuchtenDefault = TwoPMaterialLaw<Scalar, VanGenuchten, VanGenuchtenRegularization<Scalar>, TwoPEffToAbsDefaultPolicy>;

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A default configuration without regularization for using the VanGenuchten material law
 */
template<typename Scalar = double>
using VanGenuchtenNoReg = TwoPMaterialLaw<Scalar, VanGenuchten, NoRegularization<Scalar>, TwoPEffToAbsDefaultPolicy>;

} // end namespace Dumux::FluidMatrix

#endif
