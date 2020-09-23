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
 * relative permeability <-> saturation relations according to Brooks and Corey.
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_BROOKS_COREY_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_BROOKS_COREY_HH

// remove from here after release 3.3 /////////////
#include "brookscoreyparams.hh"

#include <algorithm>
#include <cmath>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
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
class [[deprecated("Use new material laws and FluidMatrix::BrooksCorey instead!")]] BrooksCorey
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    /*!
     * \brief The capillary pressure-saturation curve according to Brooks & Corey.
     *
     * The Brooks-Corey empirical capillary pressure <-> saturation
     * function is given by
     *
     *  \f$\mathrm{ p_C = p_e\overline{S}_w^{-1/\lambda}
     *  }\f$
     *
     * \param swe Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
                        and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Capillary pressure calculated by Brooks & Corey constitutive relation.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar pc(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        return params.pe()*pow(swe, -1.0/params.lambda());
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
     * \return Effective wetting phase saturation calculated as inverse of BrooksCorey constitutive relation.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar sw(const Params &params, Scalar pc)
    {
        using std::pow;
        using std::max;

        pc = max(pc, 0.0); // the equation below is undefined for negative pcs

        return pow(pc/params.pe(), -params.lambda());
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar endPointPc(const Params &params)
    { return params.pe(); }

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
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Partial derivative of \f$\mathrm{[p_c]}\f$ w.r.t. effective saturation according to Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dpc_dswe(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        return - params.pe()/params.lambda() * pow(swe, -1/params.lambda() - 1);
    }

    /*!
     * \brief The partial derivative of the effective
     *        saturation w.r.t. the capillary pressure according to Brooks & Corey.
     *
     * \param pc Capillary pressure \f$\mathrm{[p_c]}\f$ in \f$\mathrm{[Pa]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Partial derivative of effective saturation w.r.t. \f$\mathrm{[p_c]}\f$ according to Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dswe_dpc(const Params &params, Scalar pc)
    {
        using std::pow;
        using std::max;

        pc = max(pc, 0.0); // the equation below is undefined for negative pcs

        return -params.lambda()/params.pe() * pow(pc/params.pe(), - params.lambda() - 1);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by the Brooks-Corey
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Relative permeability of the wetting phase calculated as implied by Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar krw(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        return pow(swe, 2.0/params.lambda() + 3);
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase with regard to the wetting saturation of the
     *        medium implied by the Brooks-Corey parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the relative permeability of the wetting phase w.r.t. effective wetting phase
     *                  saturation calculated as implied by Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dkrw_dswe(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        return (2.0/params.lambda() + 3)*pow(swe, 2.0/params.lambda() + 2);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Relative permeability of the non-wetting phase calculated as implied by Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar krn(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar exponent = 2.0/params.lambda() + 1;
        const Scalar tmp = 1.0 - swe;
        return tmp*tmp*(1.0 - pow(swe, exponent));
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the relative permeability of the non-wetting phase w.r.t. effective wetting phase
     *                  saturation calculated as implied by Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dkrn_dswe(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const auto lambdaInv = 1.0/params.lambda();
        const auto swePow = pow(swe, 2*lambdaInv);
        return 2*(swe - 1.0)*(1.0 + (0.5 + lambdaInv)*swePow - (1.5 + lambdaInv)*swePow*swe);

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
class BrooksCorey
{
public:
    /*!
     * \brief The parameter type
     * \tparam Scalar The scalar type
     * \note The Brooks Corey laws are parameterized with two parameters: \f$\mathrm{p_e, \lambda}\f$,
     *       the entry pressure in \f$\mathrm{[Pa]}\f$] and a dimensionless shape parameter, respectively.
     */
    template<class Scalar>
    struct Params
    {
        Scalar pe, lambda;
    };

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    template<class Scalar = double>
    static Params<Scalar> makeParams(const std::string& paramGroup)
    {
        const auto pe = getParamFromGroup<Scalar>(paramGroup, "BcEntryPressure");
        const auto lambda = getParamFromGroup<Scalar>(paramGroup, "BcLambda");
        return {pe, lambda};
    }

    /*!
     * \brief The capillary pressure-saturation curve according to Brooks & Corey.
     *
     * The Brooks-Corey empirical capillary pressure <-> saturation
     * function is given by
     *
     *  \f$\mathrm{ p_C = p_e\overline{S}_w^{-1/\lambda}
     *  }\f$
     *
     * \param swe Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
                        and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Capillary pressure calculated by Brooks & Corey constitutive relation.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar pc(Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        return params.pe*pow(swe, -1.0/params.lambda);
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
     * \return Effective wetting phase saturation calculated as inverse of BrooksCorey constitutive relation.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar swe(Scalar pc, const Params<Scalar>& params)
    {
        using std::pow;
        using std::max;

        pc = max(pc, 0.0); // the equation below is undefined for negative pcs

        return pow(pc/params.pe, -params.lambda);
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    template<class Scalar>
    static Scalar endPointPc(const Params<Scalar>& params)
    { return params.pe; }

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
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Partial derivative of \f$\mathrm{[p_c]}\f$ w.r.t. effective saturation according to Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar dpc_dswe(Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        return -params.pe/params.lambda * pow(swe, -1.0/params.lambda - 1.0);
    }

    /*!
     * \brief The partial derivative of the effective
     *        saturation w.r.t. the capillary pressure according to Brooks & Corey.
     *
     * \param pc Capillary pressure \f$\mathrm{[p_c]}\f$ in \f$\mathrm{[Pa]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Partial derivative of effective saturation w.r.t. \f$\mathrm{[p_c]}\f$ according to Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar dswe_dpc(Scalar pc, const Params<Scalar>& params)
    {
        using std::pow;
        using std::max;

        pc = max(pc, 0.0); // the equation below is undefined for negative pcs

        return -params.lambda/params.pe * pow(pc/params.pe, - params.lambda - 1.0);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by the Brooks-Corey
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Relative permeability of the wetting phase calculated as implied by Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar krw(Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        return pow(swe, 2.0/params.lambda + 3.0);
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase with regard to the wetting saturation of the
     *        medium implied by the Brooks-Corey parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the relative permeability of the wetting phase w.r.t. effective wetting phase
     *                  saturation calculated as implied by Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar dkrw_dswe(Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        return (2.0/params.lambda + 3.0)*pow(swe, 2.0/params.lambda + 2.0);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Relative permeability of the non-wetting phase calculated as implied by Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar krn(Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar exponent = 2.0/params.lambda + 1.0;
        const Scalar sne = 1.0 - swe;
        return sne*sne*(1.0 - pow(swe, exponent));
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the relative permeability of the non-wetting phase w.r.t. effective wetting phase
     *                  saturation calculated as implied by Brooks & Corey.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    template<class Scalar>
    static Scalar dkrn_dswe(Scalar swe, const Params<Scalar>& params)
    {
        using std::pow;
        using std::clamp;

        swe = clamp(swe, 0.0, 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const auto lambdaInv = 1.0/params.lambda;
        const auto swePow = pow(swe, 2*lambdaInv);
        return 2.0*(swe - 1.0)*(1.0 + (0.5 + lambdaInv)*swePow - (1.5 + lambdaInv)*swePow*swe);
    }
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A regularization for the BrooksCorey material law
 * \note Regularization can be turned of by setting the threshold parameters
 *       out of range (runtime) or by replacing
 *       this class by NoRegularization (compile time).
 */
template <class Scalar>
class BrooksCoreyRegularization
{
public:
    //! Regularization parameters
    template<class S>
    struct Params
    {
        S pcLowSw;
    };

    //! Initialize the spline
    template<class MaterialLaw>
    void init(const MaterialLaw* m, const std::string& paramGroup)
    {
        pcLowSw_ = getParamFromGroup<Scalar>(paramGroup, "BcPcLowSwThreshold", 0.01);
        entryPressure_ = getParamFromGroup<Scalar>(paramGroup, "BcEntryPressure");

        initPcParameters_(m, pcLowSw_);
    }

    template<class MaterialLaw, class BaseParams, class EffToAbsParams>
    void init(const MaterialLaw* m, const BaseParams& bp, const EffToAbsParams& etap, const Params<Scalar>& p)
    {
        pcLowSw_ = p.pcLowSw;
        entryPressure_ = bp.pe;

        initPcParameters_(m, pcLowSw_);
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
        if (swe <= pcLowSw_)
            return pcLowSwPcValue_ + pcDerivativeLowSw_*(swe - pcLowSw_);

        else if (swe > 1.0)
            return pcDerivativeHighSwEnd_*(swe - 1.0) + entryPressure_;

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

        else
            return {}; // no regularization
    }

    /*!
     * \brief The regularized saturation-capillary pressure curve
     */
    OptionalScalar<Scalar> swe(const Scalar pc) const
    {
        if (pc <= entryPressure_)
            return 1.0 + (pc - entryPressure_)/pcDerivativeHighSwEnd_;

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
        if (pc <= entryPressure_)
            return 1.0/pcDerivativeHighSwEnd_;

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
        if (swe <= 0.0)
            return 0.0;
        else if (swe >= 1.0)
            return 1.0;
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
        else if (swe >= 1.0)
            return 0.0;
        else
            return {}; // no regularization
    }

    /*!
     * \brief The regularized relative permeability for the non-wetting phase
     */
    OptionalScalar<Scalar> krn(const Scalar swe) const
    {
        if (swe <= 0.0)
            return 1.0;
        else if (swe >= 1.0)
            return 0.0;
        else
            return {}; // no regularization
    }

    /*!
     * \brief The regularized derivative of the relative permeability for the non-wetting phase w.r.t. saturation
     */
    OptionalScalar<Scalar> dkrn_dswe(const Scalar swe) const
    {
        if (swe <= 0.0)
            return 0.0;
        else if (swe >= 1.0)
            return 0.0;
        else
            return {}; // no regularization
    }

private:
    template<class MaterialLaw>
    void initPcParameters_(const MaterialLaw* m, const Scalar lowSw)
    {
        pcDerivativeLowSw_ = m->template dpc_dswe<false>(lowSw);
        pcDerivativeHighSwEnd_ = m->template dpc_dswe<false>(1.0);
        pcLowSwPcValue_ = m->template pc<false>(lowSw);
    }

    Scalar pcLowSw_;
    Scalar pcLowSwPcValue_;
    Scalar entryPressure_;
    Scalar pcDerivativeLowSw_;
    Scalar pcDerivativeHighSwEnd_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A default configuration for using the VanGenuchten material law
 */
template<typename Scalar = double>
using BrooksCoreyDefault = TwoPMaterialLaw<Scalar, BrooksCorey, BrooksCoreyRegularization<Scalar>, TwoPEffToAbsDefaultPolicy>;

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A default configuration without regularization for using the VanGenuchten material law
 */
template<typename Scalar = double>
using BrooksCoreyNoReg = TwoPMaterialLaw<Scalar, BrooksCorey, NoRegularization<Scalar>, TwoPEffToAbsDefaultPolicy>;

} // end namespace Dumux::FluidMatrix

#endif
