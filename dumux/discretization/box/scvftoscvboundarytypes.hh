// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
template<class BoundaryTypes, class DiscretizationMethod>
class ScvfToScvBoundaryTypes
{
public:
    ScvfToScvBoundaryTypes() = default;

    template<class Problem>
    void computeBoundaryTypes(const Problem& problem)
    {
        // only do something for box
        if (DiscretizationMethod{} == DiscretizationMethods::box)
        {
            const auto& gridGeometry = problem.gridGeometry();
            scvBoundaryTypes.resize(gridGeometry.vertexMapper().size());
            // set all equations to Neumann by default
            for (std::size_t vIdx = 0; vIdx < scvBoundaryTypes.size(); vIdx++)
                scvBoundaryTypes[vIdx].setAllNeumann();

            auto fvGeometry = localView(gridGeometry);
            for (const auto& element : elements(gridGeometry.gridView()))
            {
                // iterate over the scvfs
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
        if (DiscretizationMethod{} == DiscretizationMethods::box)
            return scvBoundaryTypes[scv.dofIndex()];
        else
            DUNE_THROW(Dune::InvalidStateException, "Only use this for the box discretization!");
    }

private:
    std::vector<BoundaryTypes> scvBoundaryTypes;
};

} // end namespace Dumux

#endif
