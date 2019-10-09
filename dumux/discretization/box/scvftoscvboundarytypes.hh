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
 * \ingroup BoxDiscretization
 * \brief Convert intersection boundary types to vertex boundary types
 */
#ifndef DUMUX_SCVF_TO_SCV_BCTYPES_HH
#define DUMUX_SCVF_TO_SCV_BCTYPES_HH

#include <vector>
#include <dune/common/exceptions.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup BoxDiscretization
 * \brief Convert intersection boundary types to vertex boundary types
 */
template<class BoundaryTypes, DiscretizationMethod discMethod>
class ScvfToScvBoundaryTypes
{
public:
    ScvfToScvBoundaryTypes() = default;

    template<class Problem>
    void computeBoundaryTypes(const Problem& problem)
    {
        // only do something for box
        if (discMethod == DiscretizationMethod::box)
        {
            const auto& gridGeometry = problem.gridGeometry();
            scvBoundaryTypes.resize(gridGeometry.vertexMapper().size());
            // set all equations to Neumann by default
            for (std::size_t vIdx = 0; vIdx < scvBoundaryTypes.size(); vIdx++)
                scvBoundaryTypes[vIdx].setAllNeumann();

            for (const auto& element : elements(gridGeometry.gridView()))
            {
                // iterate over the scvfs
                auto fvGeometry = localView(gridGeometry);
                fvGeometry.bindElement(element);

                for (const auto& scvf : scvfs(fvGeometry))
                {
                    if (!scvf.boundary())
                        continue;

                    const auto bcTypes = problem.boundaryTypes(element, scvf);
                    if (!bcTypes.hasDirichlet())
                        continue;

                    // get the inside scv belonging to this scvf and set possible Dirichlet boundary conditions
                    const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
                    for (int pvIdx = 0; pvIdx < bcTypes.size(); ++pvIdx)
                        if (bcTypes.isDirichlet(pvIdx))
                            scvBoundaryTypes[scv.dofIndex()].setDirichlet(pvIdx);
                }
            }
        }
    }

    //! get the boundary types of the scv
    template<class SubControlVolume>
    const BoundaryTypes& boundaryTypes(const SubControlVolume& scv) const
    {
        if (discMethod == DiscretizationMethod::box)
            return scvBoundaryTypes[scv.dofIndex()];
        else
            DUNE_THROW(Dune::InvalidStateException, "Only use this for the box discretization!");
    }

private:
    std::vector<BoundaryTypes> scvBoundaryTypes;
};

} // end namespace Dumux

#endif
