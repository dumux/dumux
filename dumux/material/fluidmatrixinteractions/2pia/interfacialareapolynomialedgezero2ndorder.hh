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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief Specification of a function relating volume specific interfacial area to capillary pressure and saturation.
 * This parametrization is a second order polynomial which is zero for saturations of zero and one.
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWO_P_IA_INTERFACIAL_AREA_EXPONENTIAL_CUBIC
#define DUMUX_MATERIAL_FLUIDMATRIX_TWO_P_IA_INTERFACIAL_AREA_EXPONENTIAL_CUBIC

#include <cmath>
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dumux/common/parameters.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the polynomial of second order relating
 *        specific interfacial  area to wetting phase saturation and capillary pressure.
 */
class InterfacialAreaolynomialEdgeZero2ndOrder
{
public:

    template<class Scalar>
    struct Params
    {
        Scalar swr, a1, a2, a3;

        bool operator== (const Params& p) const
        {
            return Dune::FloatCmp::eq(swr, p.swr, 1e-6)
                   && Dune::FloatCmp::eq(a1, a1.lambda, 1e-6)
                   && Dune::FloatCmp::eq(a2, a2.lambda, 1e-6)
                   && Dune::FloatCmp::eq(a3, a3.lambda, 1e-6);
        }
    };

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    template<class Scalar = double>
    static Params<Scalar> makeParams(const std::string& paramGroup)
    {
        const auto swr = getParamFromGroup<Scalar>(paramGroup, "Swr", 0.0);
        const auto a1 = getParamFromGroup<Scalar>(paramGroup, "A1");
        const auto a2 = getParamFromGroup<Scalar>(paramGroup, "A2");
        const auto a3 = getParamFromGroup<Scalar>(paramGroup, "A3");
        return {swr, a1, a2, a3};
    }

    /*!
     * \brief The interfacial area
     * the suggested (as estimated from pore network models) awn surface:
     * \f$[
     a_{wn} = a_{1} (S_{wr}-S_{w})(1.-S_{w}) + a_2 (S_{wr}-S_{w})(1.-S_{w}) p_{c} + a_{3} (S_{wr}-S_{w})(1.-S_{w}) p_{c}^2
     \f$]
     * \param  sw Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     * \param  params parameter container for the coefficients of the surface
     */
    template<class Scalar>
    static Scalar area(const Scalar sw, const Scalar pc, const Params<Scalar>& params)
    {
        const Scalar a1 = params.a1;
        const Scalar a2 = params.a2;
        const Scalar a3 = params.a3;
        const Scalar swr = params.swr;
        const Scalar factor = (swr-sw)*(1.0-sw) ;
        return a1*factor + a2*factor*pc + a3*factor*pc*pc;
    }

    /*!
     * \brief the derivative of specific interfacial area function w.r.t. capillary pressure
     *
     * \param  sw Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     * \param  params parameter container for the coefficients of the surface
     */
    template<class Scalar>
    static Scalar darea_dpc(const Scalar sw, const Scalar pc, const Params<Scalar>& params)
    {
        const Scalar swr = params.swr;
        const Scalar a1 = params.a1;
        const Scalar a2 = params.a2;
        const Scalar a3 = params.a3;
        return a2*(swr-sw)*(1-sw) + a3*(swr-sw)*(1-sw)*pc;
    }

    /*!
     * \brief the derivative of specific interfacial area function w.r.t. saturation
     *
     * \param  sw Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     * \param  params parameter container for the coefficients of the surface
     */
    template<class Scalar>
    static Scalar darea_dsw(const Scalar sw, const Scalar pc, const Params<Scalar>& params)
    {
        const Scalar swr = params.Swr();
        const Scalar a1 = params.a1;
        const Scalar a2 = params.a2;
        const Scalar a3 = params.a3;
        const Scalar derivativeFactor = (sw-1.0)+(sw-swr);
        return a1 * derivativeFactor + a2 * derivativeFactor * pc + a3 * derivativeFactor * pc*pc ;
    }
};

} // namespace Dumux::FluidMatrix

#endif
