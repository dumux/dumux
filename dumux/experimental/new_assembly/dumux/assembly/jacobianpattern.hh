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
 * \ingroup Assembly
 * \brief Free functions to get the jacobian patterns of discretization schemes.
 */
#ifndef DUMUX_ASSEMBLY_JACOBIAN_PATTERN_HH
#define DUMUX_ASSEMBLY_JACOBIAN_PATTERN_HH

#include <dumux/discretization/method.hh>
#include <dumux/experimental/new_assembly/dumux/assembly/matrixpattern.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief Get the matrix pattern for a grid geometry using TPFA discretization.
 */
template<typename GridGeometry>
MatrixPattern getJacobianPattern(const GridGeometry& gridGeometry,
                                 const DiscretizationMethods::CCTpfa&)
{
    MatrixPattern pattern(gridGeometry.numDofs(), gridGeometry.numDofs());
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        const auto fvGeometry = localView(gridGeometry).bind(element);
        for (const auto& scv : scvs(fvGeometry))
            pattern.add(scv.dofIndex(), scv.dofIndex());
        for (const auto& scvf : scvfs(fvGeometry))
            if (!scvf.boundary())
                pattern.add(scvf.insideScvIdx(), scvf.outsideScvIdx());
    }
    return pattern;
}

/*!
 * \ingroup Assembly
 * \brief Get the matrix pattern for a grid geometry.
 */
template<typename GridGeometry>
MatrixPattern getJacobianPattern(const GridGeometry& gridGeometry)
{ return getJacobianPattern(gridGeometry, GridGeometry::discMethod); }

} // namespace Dumux

#endif
