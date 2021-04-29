// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Sequential
 * \brief Class to specify the type of a boundary.
 */
#ifndef DUMUX_SEQUENTIAL_BOUNDARY_TYPES_HH
#define DUMUX_SEQUENTIAL_BOUNDARY_TYPES_HH

#include <algorithm>
#include <array>

#include <dumux/common/boundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup Sequential
 * \brief Class to specify the type of a boundary.
 */
template <int numEq>
class SequentialBoundaryTypes : public BoundaryTypes<numEq>
{
public:
    SequentialBoundaryTypes()
    : BoundaryTypes<numEq>()
    { reset(); }

    /*!
     * \brief Reset the boundary types for one equation.
     */
    void resetEq(int eqIdx)
    {
        BoundaryTypes<numEq>::resetEq(eqIdx);
        seqBoundaryInfo_[eqIdx].isOutflow = false;
    }

    /*!
     * \brief Reset the boundary types for all equations.
     *
     * After this method no equations will be disabled and neither
     * Neumann nor Dirichlet conditions will be evaluated. This
     * corresponds to a Neumann zero boundary.
     */
    void reset()
    {
        for (int eqIdx=0; eqIdx < numEq; ++eqIdx)
            resetEq(eqIdx);
    }

    /*!
     * \brief Set all boundary conditions to Neumann.
     */
    void setAllOutflow()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            setOutflow(eqIdx);
    }

    /*!
     * \brief Set a Neumann boundary condition for a single equation.
     * \param eqIdx The index of the equation on which the outflow
     *              condition applies.
     */
    void setOutflow(int eqIdx)
    {
        resetEq(eqIdx);
        this->boundaryInfo_[eqIdx].visited = true;
        seqBoundaryInfo_[eqIdx].isOutflow = true;
    }

    /*!
     * \brief Returns true if an equation is used to specify an
     *        outflow condition.
     *
     * \param eqIdx The index of the equation
     */
    bool isOutflow(unsigned eqIdx) const
    { return seqBoundaryInfo_[eqIdx].isOutflow; }

    /*!
     * \brief Returns true if some equation is used to specify an
     *        outflow condition.
     */
    bool hasOutflow() const
    {
        return std::any_of(seqBoundaryInfo_.begin(),
                           seqBoundaryInfo_.end(),
                           [](const BoundaryInfo& b){ return b.isOutflow; }
                           );
    }

protected:
    //! use bitfields to minimize the size
    struct BoundaryInfo {
        bool isOutflow : 1;
    };

    std::array<BoundaryInfo, numEq> seqBoundaryInfo_;
};

} // end namespace Dumux

#endif
