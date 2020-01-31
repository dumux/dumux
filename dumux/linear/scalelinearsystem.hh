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
 * \ingroup Linear
 * \brief Free function to scale a linear system
 */
#ifndef DUMUX_LINEAR_SCALELINEARSYSTEM_HH
#define DUMUX_LINEAR_SCALELINEARSYSTEM_HH

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Scale the linear system by the inverse of
 * its (block-)diagonal entries.
 *
 * \param matrix the matrix to scale
 * \param rhs the right hand side vector to scale
 */
template <class Matrix, class Vector>
void scaleLinearSystem(Matrix& matrix, Vector& rhs)
{
    for (auto rowIt = matrix.begin(); rowIt != matrix.end(); ++rowIt)
    {
        auto rowIdx = rowIt.index();
        auto diagonal = matrix[rowIdx][rowIdx];
        diagonal.invert();

        const auto b = rhs[rowIdx];
        diagonal.mv(b, rhs[rowIdx]);

        for (auto& col : *rowIt)
            col.leftmultiply(diagonal);
    }
}

} // end namespace Dumux
