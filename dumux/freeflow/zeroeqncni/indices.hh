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
 * \brief Defines the indices required for the non-isothermal
 *        compositional ZeroEq box model.
 */
#ifndef DUMUX_ZEROEQNCNI_INDICES_HH
#define DUMUX_ZEROEQNCNI_INDICES_HH

#include <dumux/freeflow/stokesncni/indices.hh>
#include <dumux/freeflow/zeroeq/indices.hh>
#include <dumux/freeflow/zeroeqnc/indices.hh>

namespace Dumux
{

/*!
 * \ingroup BoxZeroEqncniModel
 * \brief The indices for the eddy conductivity model.
 */
struct EddyConductivityIndices
{
    // Eddy Conductivity Model Indices
    static const int noEddyConductivityModel = 0;
    static const int reynoldsAnalogy = 1;
    static const int modifiedVanDriest = 2;
    static const int deissler = 3;
    static const int meier = 4;
};

/*!
 * \ingroup BoxZeroEqncniModel
 * \ingroup ImplicitIndices
 * \brief The common indices for the non-isothermal compositional ZeroEq box model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset = 0>
struct ZeroEqncniCommonIndices : public StokesncniCommonIndices<TypeTag, PVOffset>
{
    static const int scvDataPrecision = 5;
    static const int scvDataWidth = scvDataPrecision + 10;
};

} // end namespace

#endif // DUMUX_ZEROEQNCNI_INDICES_HH
