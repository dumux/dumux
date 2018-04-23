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
 * \ingroup FreeflowNCModel
 * \copydoc Dumux::FreeflowNCIndices
 */
#ifndef DUMUX_FREEFLOW_NC_INDICES_HH
#define DUMUX_FREEFLOW_NC_INDICES_HH

#include <dumux/freeflow/navierstokes/indices.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowNCModel
 * \brief The common indices for the isothermal multi-component free-flow model.
 */
template <int dimension, int numEquations,
          int phaseIdx, int theReplaceCompEqIdx>
struct FreeflowNCIndices : public NavierStokesIndices<dimension>
{
public:
    //! The index of the fluid phase in the fluid system
    static constexpr int fluidSystemPhaseIdx = phaseIdx;

    //! The index of the main component
    static constexpr int mainCompIdx = fluidSystemPhaseIdx;

    //! The index of the component whose mass balance will be replaced by the total one
    static constexpr int replaceCompEqIdx = theReplaceCompEqIdx;
};

} // end namespace Dumux

#endif
