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
#ifndef LEVERETT_LAW_HH
#define LEVERETT_LAW_HH

#include "leverettlawparams.hh"

#include <dumux/common/spline.hh>

#include <algorithm>

#include <math.h>
#include <assert.h>

namespace Dumux
{
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the capillary pressure <-> saturation
 *        relation for the heatpipe problem.
 *
 * This class bundles the "raw" curves as static members and doesn't concern itself
 * converting absolute to effective saturations and vince versa.
 */
template <class BaseLaw, class ParamsT = LeverettLawParams<typename BaseLaw::Params> >
class LeverettLaw : public BaseLaw
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
    static Scalar pc(const Params &params, Scalar swe)
    {
        return BaseLaw::pc(params, swe) * params.leverettFactor();

    }



private:

};

}

#endif
