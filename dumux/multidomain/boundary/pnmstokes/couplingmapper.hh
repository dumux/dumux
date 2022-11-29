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
 * \ingroup BoubdaryCoupling
 * \ingroup BoxModel
 * \copydoc Dumux::PNMStokesCouplingManager
 */

#ifndef DUMUX_PNM_STOKES_COUPLINGMAPPER_HH
#define DUMUX_PNM_STOKES_COUPLINGMAPPER_HH

#include <memory>
#include <vector>
#include <type_traits>

#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dumux/multidomain/boundary/pnmstokes/geometry.hh>
#include <dumux/geometry/intersectingentities.hh>

namespace Dumux {

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionBoundary
 * \brief Coupling manager for low-dimensional domains coupled at the boundary to the bulk
 *        domain.
 */
class PNMStokesCouplingMapper
{

public:

    //! retruns the relevant coupling maps and stencils
    template<class CouplingManager>
    static auto computeCouplingMapsAndStencils(const CouplingManager& couplingManager)
    {
        using std::abs;
        struct PNMStokesMappingData
        {
            using Stencils = typename CouplingManager::CouplingStencils;
            Stencils lowDimToBulkCellCenterStencils;
            Stencils lowDimToBulkFaceStencils;
            Stencils bulkCellCenterToLowDimStencils;
            Stencils bulkFaceToLowDimStencils;

            std::vector<bool> isCoupledLowDimDof;
            std::vector<bool> isCoupledBulkFaceDof;
            std::vector<bool> isCoupledBulkFrontalFaceDof;

            std::unordered_map<std::size_t, std::vector<std::size_t>> lowDimElementToBulkElementsMap;
            std::unordered_map<std::size_t, std::size_t> bulkElementToLowDimElementMap;
        } data;

        const auto& bulkGridGeometry = couplingManager.problem(CouplingManager::bulkIdx).gridGeometry();
        const auto& lowDimGridGeometry = couplingManager.problem(CouplingManager::lowDimIdx).gridGeometry();

        using BulkGridView = std::decay_t<decltype(bulkGridGeometry.gridView())>;
        using GlobalPosition = typename BulkGridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
        using Scalar = typename BulkGridView::ctype;

        data.isCoupledLowDimDof.resize(lowDimGridGeometry.numDofs(), false);
        data.isCoupledBulkFaceDof.resize(bulkGridGeometry.numFaceDofs(), false);
        data.isCoupledBulkFrontalFaceDof.resize(bulkGridGeometry.numFaceDofs(), false);

        auto lowDimFvGeometry = localView(lowDimGridGeometry);
        auto bulkFvGeometry = localView(bulkGridGeometry);

        using ThroatInterfaceGeometryType = Dune::AxisAlignedCubeGeometry<Scalar,
                                                                          BulkGridView::dimension-1,
                                                                          BulkGridView::dimensionworld>;

        // iterate over the lowDim elements
        for (const auto& lowDimElement : elements(lowDimGridGeometry.gridView()))
        {
            // check if the lowDim vertices intersect with the bulk grid
            lowDimFvGeometry.bindElement(lowDimElement);
            for (const auto& lowDimScv : scvs(lowDimFvGeometry))
            {
                // skip the dof if it is not on the boundary
                if (!lowDimGridGeometry.dofOnBoundary(lowDimScv.dofIndex()))
                    continue;

                // get the intersection bulk element
                const auto lowDimPos = lowDimScv.dofPosition();
                const auto lowDimDofIdx = lowDimScv.dofIndex();

                // check for intersections, skip if no intersection was found
                if (intersectingEntities(lowDimPos, bulkGridGeometry.boundingBoxTree()).empty())
                    continue;
                else
                    data.isCoupledLowDimDof[lowDimDofIdx] = true;

                const auto lowDimElementIdx = lowDimGridGeometry.elementMapper().index(lowDimElement);
                const auto& otherLowDimScv = lowDimFvGeometry.scv(1 - lowDimScv.indexInElement());
                const auto otherLowDimScvDofIdx = otherLowDimScv.dofIndex();

                static const auto couplingPlaneNormal = getParamFromGroup<GlobalPosition>(couplingManager.problem(CouplingManager::bulkIdx).paramGroup(),
                                                                                          "Grid.CouplingPlaneNormal",
                                                                                          [](){ GlobalPosition tmp(0.0); tmp[tmp.size()-1] = 1.0; return tmp; }());

                // Check whether to use the pore body radius or the projected throat radius for coupling.
                // If needed, determine the throat radius, which might depend on the throat's angle of orientation.
                static const bool coupleOverPoreRadius = getParamFromGroup<bool>(couplingManager.problem(CouplingManager::bulkIdx).paramGroup(), "Grid.CoupleOverPoreRadius", false);
                const Scalar couplingThroatRadius = coupleOverPoreRadius ? lowDimGridGeometry.poreInscribedRadius(lowDimDofIdx)
                                                  : projectedThroatRadius(lowDimGridGeometry.throatInscribedRadius(lowDimElementIdx),
                                                                          lowDimElement, couplingPlaneNormal);

                using std::abs;
                const int couplingPlaneDirectionIdx = std::find_if(couplingPlaneNormal.begin(), couplingPlaneNormal.end(),
                                                                   [eps = couplingThroatRadius*1e-8](const auto& x) { return abs(x) > eps; } ) - couplingPlaneNormal.begin();

                auto lowerLeft = lowDimPos - GlobalPosition(couplingThroatRadius);
                lowerLeft[couplingPlaneDirectionIdx] = lowDimPos[couplingPlaneDirectionIdx];
                auto upperRight = lowDimPos + GlobalPosition(couplingThroatRadius);
                upperRight[couplingPlaneDirectionIdx] = lowDimPos[couplingPlaneDirectionIdx];

                auto axes = std::move(std::bitset<BulkGridView::dimensionworld>{}.set());
                axes.set(couplingPlaneDirectionIdx, false);

                ThroatInterfaceGeometryType throatGeometry(lowerLeft, upperRight, axes);
                const auto interfaceIntersections = intersectingEntities(std::move(throatGeometry), bulkGridGeometry.boundingBoxTree());

                for (const auto intersection : interfaceIntersections)
                {
                    const auto bulkElemIdx = intersection.second();
                    data.lowDimElementToBulkElementsMap[lowDimElementIdx].push_back(bulkElemIdx);

                    if (!data.bulkElementToLowDimElementMap.count(bulkElemIdx))
                            data.bulkElementToLowDimElementMap[bulkElemIdx] = lowDimElementIdx;

                    data.lowDimToBulkCellCenterStencils[lowDimElementIdx].push_back(bulkElemIdx);
                    data.bulkCellCenterToLowDimStencils[bulkElemIdx].push_back(lowDimDofIdx);

                    const auto& bulkElement = bulkGridGeometry.boundingBoxTree().entitySet().entity(bulkElemIdx);
                    bulkFvGeometry.bindElement(bulkElement);

                    const auto coupledFaceDofIndices = coupledFaces_(bulkFvGeometry, lowDimPos, couplingThroatRadius, couplingPlaneDirectionIdx);

                    data.lowDimToBulkFaceStencils[lowDimElementIdx].push_back(coupledFaceDofIndices.coupledFrontalFace);
                    data.bulkFaceToLowDimStencils[coupledFaceDofIndices.coupledFrontalFace].push_back(lowDimDofIdx);
                    data.bulkFaceToLowDimStencils[coupledFaceDofIndices.coupledFrontalFace].push_back(otherLowDimScvDofIdx);

                    data.isCoupledBulkFaceDof[coupledFaceDofIndices.coupledFrontalFace] = true;
                    data.isCoupledBulkFrontalFaceDof[coupledFaceDofIndices.coupledFrontalFace] = true;

                        // treat the coupled normal faces
                    for (const auto bulkFaceDofIdx : coupledFaceDofIndices.coupledNormalFaces)
                    {
                        data.bulkFaceToLowDimStencils[bulkFaceDofIdx].push_back(lowDimDofIdx);
                        data.bulkFaceToLowDimStencils[bulkFaceDofIdx].push_back(otherLowDimScvDofIdx);
                        data.isCoupledBulkFaceDof[bulkFaceDofIdx] = true;
                    }
                }
            }
        }

        return data;
    }


private:

    //! get the bulk faces that are coupled to the lowDim dofs
    template<class FVElementGeometry, class GlobalPosition, class Scalar>
    static auto coupledFaces_(const FVElementGeometry& fvGeometry,
                              const GlobalPosition& lowDimPos,
                              const Scalar couplingThroatRadius,
                              const int couplingPlaneDirectionIdx)
    {
        using std::abs;
        static constexpr auto dimWorld = GlobalPosition::dimension;

        struct Result
        {
            std::vector<std::size_t> coupledNormalFaces;
            std::size_t coupledFrontalFace;
        } result;

        result.coupledNormalFaces.reserve(fvGeometry.numScvf());

        for (const auto& scvf : scvfs(fvGeometry))
        {
            const Scalar eps = scvf.area()*1e-8;
            assert(eps < couplingThroatRadius);

            if (scvf.directionIndex() == couplingPlaneDirectionIdx) // the bulk faces that lie within the coupling interface
            {
                if (scvf.center()[couplingPlaneDirectionIdx] - lowDimPos[couplingPlaneDirectionIdx] < eps)
                    result.coupledFrontalFace = scvf.dofIndex();
            }
            else // the bulk faces perpendicular to the coupling interface
            {
                bool isCoupledFace = false;

                for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                {
                    if (dimIdx == couplingPlaneDirectionIdx || scvf.boundary())
                        continue;

                    isCoupledFace = abs(scvf.center()[dimIdx] - lowDimPos[dimIdx]) < couplingThroatRadius + eps;
                }

                if (isCoupledFace)
                    result.coupledNormalFaces.push_back(scvf.dofIndex());
            }
        }
        return result;
    }
};

} // end namespace Dumux

#endif
