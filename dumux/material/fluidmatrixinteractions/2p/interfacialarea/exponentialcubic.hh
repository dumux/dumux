// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief Specification of a function relating volume specific interfacial area to capillary pressure and saturation.
 *        This function is of third order in pc.
 * \note  It is used for calculating the interfacial area between the nonwetting and solid phase
 *        by Nuske 2014 (https://elib.uni-stuttgart.de/handle/11682/614, page 62) \cite nuske2014.
 *
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWO_P_INTERFACIAL_AREA_EXPONENTIAL_CUBIC
#define DUMUX_MATERIAL_FLUIDMATRIX_TWO_P_INTERFACIAL_AREA_EXPONENTIAL_CUBIC

#include <cmath>
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dumux/common/parameters.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of a exponential function relating
 *        specific interfacial  area to wetting phase saturation and capillary pressure.
 */
class InterfacialAreaExponentialCubic
{
public:

    template<class Scalar>
    struct Params
    {
        Params(Scalar a1 = 0, Scalar a2 = 0, Scalar a3 = 0)
        : a1_(a1), a2_(a2), a3_(a3) {}

        Scalar a1() const { return a1_; }
        void setA1(Scalar a1) { a1_ = a1; }

        Scalar a2() const { return a2_; }
        void setA2(Scalar a2) { a2_ = a2; }

        Scalar a3() const { return a3_; }
        void setA3(Scalar a3) { a3_ = a3; }

        bool operator== (const Params& p) const
        {
            return Dune::FloatCmp::eq(a1(), p.a1(), 1e-6)
                   && Dune::FloatCmp::eq(a2(), p.a2(), 1e-6)
                   && Dune::FloatCmp::eq(a3(), p.a3(), 1e-6);
        }

    private:
        Scalar a1_, a2_, a3_;
    };

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    template<class Scalar = double>
    static Params<Scalar> makeParams(const std::string& paramGroup)
    {
        const auto a1 = getParamFromGroup<Scalar>(paramGroup, "A1");
        const auto a2 = getParamFromGroup<Scalar>(paramGroup, "A2");
        const auto a3 = getParamFromGroup<Scalar>(paramGroup, "A3");
        return {a1, a2, a3};
    }

    /*!
     * \brief The interfacial area
     *
     * the suggested (as estimated from pore network models) interfacial area between the nonwetting and solid phase:
     * \f$\mathrm{
        a_{ns} = a_1  e^{a_2 * S_w }  + a_3 * p_c^3 ;
     }\f$
     * \param  swe Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     * \param  params parameter container for the coefficients of the surface
     */
    template<class Scalar>
    static Scalar area(const Scalar swe, const Scalar pc, const Params<Scalar>& params)
    {
        // TODO think about awn surface for relative saturation
        const Scalar a1 = params.a1();
        const Scalar a2 = params.a2();
        const Scalar a3 = params.a3();

        using std::exp;
        return a1 * exp( a2 * swe) + a3 * pc * pc * pc ;
    }
};

} // namespace Dumux::FluidMatrix

#endif
