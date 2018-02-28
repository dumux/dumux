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
 * \ingroup NavierStokesNCModel
 * \copydoc Dumux::NavierStokesNCIndices
 */
#ifndef DUMUX_STAGGERED_NAVIERSTOKES_NC_INDICES_HH
#define DUMUX_STAGGERED_NAVIERSTOKES_NC_INDICES_HH

#include <dumux/freeflow/navierstokes/indices.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

// \{
/*!
 * \ingroup NavierStokesNCModel
 * \brief The common indices for the isothermal multi-component Navier-Stokes model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int dimension, int numEquations,
          int thePhaseIdx, int theReplaceCompEqIdx,
          int PVOffset = 0>
struct NavierStokesNCIndices : public NavierStokesIndices<dimension, numEquations, PVOffset>
{
private:
    using ParentType = NavierStokesIndices<dimension, numEquations, PVOffset>;

public:

    static constexpr int phaseIdx = thePhaseIdx; //!< The phase index
    static constexpr int mainCompIdx = phaseIdx; //!< The index of the main component

    //! The index of the component whose mass balance will be replaced by the total one
    static constexpr int replaceCompEqIdx = theReplaceCompEqIdx;
};

// \}
} // end namespace

#endif
