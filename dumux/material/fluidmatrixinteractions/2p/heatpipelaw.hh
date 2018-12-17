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
 * \brief Implementation of the capillary pressure <-> saturation relation
 *        for the heatpipe problem.
 */
#ifndef HEATPIPELAW_HH
#define HEATPIPELAW_HH

#include "heatpipelawparams.hh"

#include <dumux/common/spline.hh>

#include <algorithm>

#include <math.h>
#include <assert.h>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the capillary pressure <-> saturation
 *        relation for the heatpipe problem.
 *
 * This class bundles the "raw" curves as static members and doesn't concern itself
 * converting absolute to effective saturations and vince versa.
 */
template <class ScalarT, class ParamsT = HeatPipeLawParams<ScalarT> >
class HeatPipeLaw
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     * \param params Array of parameters asd
     * \param Sw Effective saturation of of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     */
    static Scalar pc(const Params &params, Scalar Sw)
    {
        Scalar Sn = 1 - Sw;
        Scalar p0Gamma = params.p0()*params.gamma();

        // regularization
        if (Sn >= 1.0) {
            Scalar y = p0Gamma*(  (1.263*1.0 -   2.120)*1.0 + 1.417)*1.0;
            Scalar m = p0Gamma*((3*1.263*1.0 - 2*2.120)*1.0 + 1.417);
            return (Sn - 1)*m + y;
        }
        else if (Sn <= 0.0) {
            Scalar y = 0.0;
            Scalar m = p0Gamma*1.417;
            return Sn*m + y;
        }

        return p0Gamma*((1.263*Sn - 2.120)*Sn + 1.417) * Sn;
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * \return The effective saturaion of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params Array of parameters
     * \param pC  capillary pressure \f$\mathrm{[p_C]}\f$ in \f$\mathrm{[Pa]}\f$.
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        DUNE_THROW(Dune::NotImplemented, "HeatPipeLaw::Sw");
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     * \param params Array of parameters
     * \param Sw Effective saturation of of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     */
    static Scalar dpC_dSw(const Params &params, Scalar Sw)
    {
        Scalar Sn = 1 - Sw;
        Scalar p0Gamma = params.p0()*params.gamma();
        if (Sn > 1.0)
            Sn = 1.0;
        else if (Sn <= 0.0) {
            Scalar m = -p0Gamma*1.417;
            return m;
        }

        Scalar m = - p0Gamma*((3*1.263*Sn - 2*2.120)*Sn + 1.417);
        return m;
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     * \param params Array of parameters
     * \param pC  capillary pressure \f$\mathrm{[p_C]}\f$ in \f$\mathrm{[Pa]}\f$.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        DUNE_THROW(Dune::NotImplemented, "HeatPipeLaw::dSw_dpC");
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param params Array of parameters
     * \param Sw The mobile saturation of the wetting phase.
     */
    static Scalar krw(const Params &params, Scalar Sw)
    {
        return kr_(Sw);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     * \param params Array of parameters
     * \param Sw The mobile saturation of the wetting phase.
     */
    static Scalar krn(const Params &params, Scalar Sw)
    {
        Scalar Sn = 1 - Sw;
        return kr_(Sn);
    }

private:
    static Scalar kr_(Scalar S)
    {
        const Scalar eps = 0.95;
        if (S >= 1)
            return 1;
        else if (S <= 0)
            return 0;
        else if (S > eps) {
            // regularize
            using Spline = Dumux::Spline<Scalar>;
            Spline sp(eps, 1.0, // x1, x2
                      eps*eps*eps, 1, // y1, y2
                      3*eps*eps, 0); // m1, m2
            return sp.eval(S);
        }

        return S*S*S;
    }
};

} // end namespace Dumux

#endif
