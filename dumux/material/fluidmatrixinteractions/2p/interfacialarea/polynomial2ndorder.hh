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
 * This parametrization is a second order polynomial.
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWO_P_INTERFACIAL_AREA_POLYNOMIAL_SECOND_ORDER_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_TWO_P_INTERFACIAL_AREA_POLYNOMIAL_SECOND_ORDER_HH

#include <cmath>
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dumux/common/parameters.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the polynomial of second order relating
 *        specific interfacial  area to wetting phase saturation and capillary pressure as suggested by Joekar-Niasar(2008) \cite joekar2008 .
 */
class InterfacialAreaPolynomialSecondOrder
{
public:

    template<class Scalar>
    struct Params
    {
        Params(Scalar a00 = 0, Scalar a01 = 0, Scalar a02 = 0, Scalar a10 = 0, Scalar a11 = 0, Scalar a20 = 0)
        : a00_(a00), a01_(a01), a02_(a02), a10_(a10), a11_(a11), a20_(a20) {}

        Scalar a00() const { return a00_; }
        void setA00(Scalar a00) { a00_ = a00; }

        Scalar a01() const { return a01_; }
        void setA01(Scalar a01) { a01_ = a01; }

        Scalar a02() const { return a02_; }
        void setA02(Scalar a02) { a02_ = a02; }

        Scalar a10() const { return a10_; }
        void setA10(Scalar a10) { a10_ = a10; }

        Scalar a11() const { return a11_; }
        void setA11(Scalar a11) { a11_ = a11; }

        Scalar a20() const { return a20_; }
        void setA20(Scalar a20) { a20_ = a20; }

        bool operator== (const Params& p) const
        {
            return Dune::FloatCmp::eq(a00(), p.a00(), 1e-6)
                   && Dune::FloatCmp::eq(a01(), p.a01(), 1e-6)
                   && Dune::FloatCmp::eq(a02(), p.a02(), 1e-6)
                   && Dune::FloatCmp::eq(a10(), p.a10(), 1e-6)
                   && Dune::FloatCmp::eq(a11(), p.a11(), 1e-6)
                   && Dune::FloatCmp::eq(a20(), p.a20(), 1e-6);
        }

    private:
        Scalar a00_, a01_, a02_, a10_, a11_, a20_;
    };

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    template<class Scalar = double>
    static Params<Scalar> makeParams(const std::string& paramGroup)
    {
        const auto a00 = getParamFromGroup<Scalar>(paramGroup, "A00", 0.0);
        const auto a01 = getParamFromGroup<Scalar>(paramGroup, "A01");
        const auto a02 = getParamFromGroup<Scalar>(paramGroup, "A02");
        const auto a10 = getParamFromGroup<Scalar>(paramGroup, "A10");
        const auto a11 = getParamFromGroup<Scalar>(paramGroup, "A11");
        const auto a20 = getParamFromGroup<Scalar>(paramGroup, "A20");
        return {a00, a01, a02, a10, a11, a20};
    }

    /*!
     * \brief The interfacial area
     *
     * the suggested (as estimated from pore network models) awn surface:
     * \f$[
     a_{wn} = a_{00} + a_{10}S_{w} + a_{20}S_{w}^2 + a_{11} S_{w} p_{c} +  a_{01} p_{c} + a_{02}p_{c}^2
     \f$]
     * \param  params parameter container for the coefficients of the surface
     * \param  swe Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    template<class Scalar>
    static Scalar area(const Scalar swe, const Scalar pc, const Params<Scalar>& params)
    {
        const Scalar a00 = params.a00();
        const Scalar a10 = params.a10();
        const Scalar a20 = params.a20();
        const Scalar a11 = params.a11();
        const Scalar a01 = params.a01();
        const Scalar a02 = params.a02();

        return a00 + a10 * swe + a20 * swe*swe + a11*swe*pc +  a01*pc + a02*pc*pc;
    }

    /*!
     * \brief the derivative of specific interfacial area function w.r.t. capillary pressure
     *
     * \param  params parameter container for the coefficients of the surface
     * \param  swe Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    template<class Scalar>
    static Scalar darea_dpc(const Scalar swe, const Scalar pc, const Params<Scalar>& params)
    {
        return params.a11()*swe + params.a01() + 2.0*params.a02() * pc;
    }

    /*!
     * \brief the derivative of specific interfacial area function w.r.t. saturation
     *
     * \param  params parameter container for the coefficients of the surface
     * \param  swe Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    template<class Scalar>
    static Scalar darea_dsw(const Scalar swe, const Scalar pc, const Params<Scalar>& params)
    {
        return params.a11()*pc + params.a10() + 2.0*params.a20()*swe;
    }
};

} // end namespace Dumux::FluidMatrix

#endif
