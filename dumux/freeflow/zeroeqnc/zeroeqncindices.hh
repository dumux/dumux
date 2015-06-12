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
 * \brief Defines the indices required for the compositional ZeroEq box model.
 */
#ifndef DUMUX_ZEROEQNC_INDICES_HH
#define DUMUX_ZEROEQNC_INDICES_HH

#include <dumux/freeflow/stokesnc/stokesncindices.hh>
#include <dumux/freeflow/zeroeq/zeroeqindices.hh>

namespace Dumux
{

/*!
 * \ingroup BoxZeroEqncModel
 * \brief The indices for the eddy diffusivity model.
 */
struct EddyDiffusivityIndices
{
    // Eddy Diffusivity Model Indices
    static const int noEddyDiffusivityModel = 0;
    static const int reynoldsAnalogy = 1;
    static const int modifiedVanDriest = 2;
    static const int deissler = 3;
    static const int meier = 4;
    static const int exponential = 5;
};

/*!
 * \ingroup BoxZeroEqncModel
 * \ingroup ImplicitIndices
 * \brief The common indices for the isothermal compositional ZeroEq box model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset = 0>
struct ZeroEqncCommonIndices : public StokesncCommonIndices<TypeTag, PVOffset>
{
    static const int scvDataPrecision = 5;
    static const int scvDataWidth = scvDataPrecision + 10;
};

} // end namespace

#endif // DUMUX_ZEROEQNC_INDICES_HH
