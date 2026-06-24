// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux-Lecture contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_KR_PC_HEATPIPE_HH
#define DUMUX_KR_PC_HEATPIPE_HH

#include <algorithm>
#include <cmath>

#include <dumux/common/parameters.hh>
#include <dumux/common/spline.hh>
#include <dumux/common/optionalscalar.hh>
#include <dumux/material/fluidmatrixinteractions/2p/materiallaw.hh>

namespace Dumux::FluidMatrix {
/*!
 *
 * \brief Implementation of the capillary pressure <-> saturation
 *        and relative permeability <-> saturation relation for the heatpipe
 *        problem. rel-perm is based on the model of Fatt and Klikoff,
 *        cap-press based on the function of Leverett.
 *
 */
class KrPcHeatPipe
{
public:
    template<class Scalar>
    struct Params
    {
        Params(Scalar swr, Scalar snr, Scalar p0)
        : swr_(swr)
        , snr_(snr)
        , p0_(p0)
        {}

        Scalar swr() const{ return swr_; }
        void setSwr(Scalar swr){ swr_ = swr; }

        Scalar snr() const { return snr_; }
        void setSnr(Scalar snr) { snr_ = snr; }

        Scalar p0() const { return p0_; }
        void setp0(Scalar p0) { p0_ = p0; }

        bool operator== (const Params& p) const
        {
            return Dune::FloatCmp::eq(snr(), p.snr(), 1e-6)
                   && Dune::FloatCmp::eq(snr(), p.snr(), 1e-6)
                   && Dune::FloatCmp::eq(p0(), p.p0(), 1e-6);
        }

    private:
        Scalar swr_, snr_, p0_;
    };

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    template<class Scalar = double>
    static Params<Scalar> makeParams(const std::string& paramGroup)
    {
        const auto swr = getParamFromGroup<Scalar>(paramGroup, "Swr");
        const auto snr = getParamFromGroup<Scalar>(paramGroup, "Snr");
        const auto p0 = getParamFromGroup<Scalar>(paramGroup, "P0");
        return {swr, snr, p0};
    }

    /*!
     * \brief The capillary pressure-saturation curve according to Leverett.
     *
     */
    template<class Scalar>
    static Scalar pc(Scalar sw, const Params<Scalar>& params)
    {
        if(sw<0.) sw=0.;
         /* effective values */
        Scalar Swe = (sw-params.swr())/(1-params.snr()-params.swr());
        if(Swe<0) Swe=0.;
        if(Swe>1.0) Swe=1.0;
        Scalar f = 1.417*(1-Swe)-2.120*std::pow((1-Swe),2)+1.263*std::pow((1-Swe),3);
        Scalar sigma=0.0588;
        Scalar value= params.p0()*sigma*f;
        /* regularization for small saturations */
        if(sw<params.swr())
        { Scalar a = -0.966*params.p0()*sigma/(4.*std::pow(params.swr(),3));
          Scalar c = 0.56*params.p0()*sigma - a*std::pow(params.swr(),4);
          value = a*std::pow(sw,4)+c;
        }

        Scalar MODIFY_ME_CAPILLARY_PRESSURE = getParam<Scalar>("Problem.MultiplierPC");

        return(value * MODIFY_ME_CAPILLARY_PRESSURE);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium according to the Fatt-Klikoff
     *        parameterization.
     */
    template<class Scalar>
    static Scalar krw(Scalar sw, const Params<Scalar>& params)
    {
         /*** effective saturation ***/
         Scalar Se = (sw-params.swr())/(1-params.snr()-params.swr());
         /* effective Saturation Se has to be between 0 and 1! */
         if (Se>1.)  Se=1.;
         if (Se<0) Se=0;
         /* compute and return value */
         return std::pow(Se,3);
    };

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium according to the Fatt-Klikoff
     *        parameterization.
     */
    template<class Scalar>
    static Scalar krn(Scalar sw, const Params<Scalar>& params)
    {
         /*** effective saturation ***/
         Scalar Se = (sw-params.swr())/(1-params.snr()-params.swr());
         /* effective Saturation Se has to be between 0 and 1! */
         if (Se>1.)  Se=1.;
         if (Se<0) Se=0;
         /* compute and return value */
         return std::pow(1.-Se,3);
    };
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A default configuration for using the VanGenuchten material law
 */
template<typename Scalar = double>
using KrPcHeatPipeDefault = TwoPMaterialLaw<Scalar, KrPcHeatPipe, NoRegularization, TwoPEffToAbsDefaultPolicy>;

}

#endif
