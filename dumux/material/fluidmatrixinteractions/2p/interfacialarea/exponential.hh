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
 * This function is exponential.
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWO_P_INTERFACIAL_AREA_EXPONENTIAL
#define DUMUX_MATERIAL_FLUIDMATRIX_TWO_P_INTERFACIAL_AREA_EXPONENTIAL

#include <cmath>
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dumux/common/parameters.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the exponential function relating
 *        specific interfacial area to wetting phase saturation and
 *        capillary pressure as suggested by Nuske(2009) (Diploma thesis) \cite nuske2009 .
 */
class InterfacialAreaExponential
{
public:

    template<class Scalar>
    struct Params
    {
        Params(Scalar swr = 0, Scalar a1 = 0, Scalar a2 = 0, Scalar a3 = 0)
        : swr_(swr), a1_(a1), a2_(a2), a3_(a3) {}

        Scalar swr() const { return swr_; }
        void setSwr(Scalar swr) { swr_ = swr; }

        Scalar a1() const { return a1_; }
        void setA1(Scalar a1) { a1_ = a1; }

        Scalar a2() const { return a2_; }
        void setA2(Scalar a2) { a2_ = a2; }

        Scalar a3() const { return a3_; }
        void setA3(Scalar a3) { a3_ = a3; }

        bool operator== (const Params& p) const
        {
            return Dune::FloatCmp::eq(swr(), p.swr(), 1e-6)
                   && Dune::FloatCmp::eq(a1(), p.a1(), 1e-6)
                   && Dune::FloatCmp::eq(a2(), p.a2(), 1e-6)
                   && Dune::FloatCmp::eq(a3(), p.a3(), 1e-6);
        }

    private:
        Scalar swr_, a1_, a2_, a3_;
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
     *
     * the suggested (as estimated from pore network models) awn surface:
     * \f$\mathrm{
        a_{wn} = a_1 * (S_{wr}-S_w) .* (1-S_w) + a_2 * (S_{wr}-S_w) * (1-S_w) * \exp( a_3 * p_c) ;
     }\f$
     * \param  swe Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     * \param  params parameter container for the coefficients of the surface
     */
    template<class Scalar>
    static Scalar area(const Scalar swe, const Scalar pc, const Params<Scalar>& params)
    {
        const Scalar a1 = params.a1();
        const Scalar a2 = params.a2();
        const Scalar a3 = params.a3();
        const Scalar swr = params.swr();
        using std::exp;
        return a1 * (swr-swe) * (1-swe) + a2 * (swr-swe) * (1-swe) * exp( a3 * pc) ;
    }

    /*! \brief the derivative of specific interfacial area function w.r.t. capillary pressure
     *
     * \param  params parameter container for the coefficients of the surface
     * \param  swe Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    template<class Scalar>
    static Scalar darea_dpc(const Scalar swe, const Scalar pc, const Params<Scalar>& params)
    {
        const Scalar a2 = params.a2();
        const Scalar a3 = params.a3();
        const Scalar swr = params.swr();
        using std::exp;
        return a2 * a3 * (swr-swe) * (1-swe) * exp(a3*pc);
    }

    /*! \brief the derivative of specific interfacial area function w.r.t. saturation
     *
     * \param  params parameter container for the coefficients of the surface
     * \param  swe Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    template<class Scalar>
    static Scalar darea_dsw(const Scalar swe, const Scalar pc, const Params<Scalar>& params)
    {
        Scalar a1 = params.a1();
        Scalar a2 = params.a2();
        Scalar a3 = params.a3();
        Scalar swr = params.swr();
        using std::exp;
        return - a1 *(1+swr-2*swe) - a2 * exp(a3*pc) * (1+swr-2*swe);
    }

};

} // namespace Dumux::FluidMatrix

#endif
