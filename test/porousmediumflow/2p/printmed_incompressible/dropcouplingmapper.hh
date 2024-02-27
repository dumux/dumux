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
 * \ingroup BoubdaryCoupling
 * \ingroup BoxModel
 * \copydoc Dumux::PNMStokesCouplingManager
 */

#ifndef DUMUX_DROP_COUPLINGMAPPER_HH
#define DUMUX_DROP_COUPLINGMAPPER_HH

#include <memory>
#include <vector>
#include <type_traits>

#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dumux/mixeddimension/boundary/pnmstokes/geometry.hh>

#include "dropintersection.hh"

namespace Dumux {

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionBoundary
 * \brief Coupling manager for drop coupled at the boundary to the bulk and low-dim
 *        domain.
 */

class DropCouplingMapper
{
    public:

    //! retruns the relevant coupling maps and stencils
    template<class CouplingManager>
    static auto computeCouplingMapsAndStencils(const CouplingManager& couplingManager)
    {

        using std::abs;
        const auto& bulkGridGeometry = couplingManager.problem(CouplingManager::bulkIdx).gridGeometry();
        const auto& dropGridGeometry = couplingManager.problem(CouplingManager::dropIdx).gridGeometry();

        auto dropFvGeometry = localView(dropGridGeometry);
        auto bulkFvGeometry = localView(bulkGridGeometry);

        using BulkGridView = std::decay_t<decltype(bulkGridGeometry.gridView())>;
        using GlobalPosition = typename BulkGridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
        using Scalar = typename BulkGridView::ctype;

        struct DropStokesMappingData
        {
            using Stencils = typename CouplingManager::CouplingStencils;
            Stencils dropToBulkCellCenterStencils;
            Stencils dropToBulkFaceStencils;

            Stencils bulkCellCenterToDropStencils;
            Stencils bulkFaceToDropStencils;

            std::vector<bool> isCoupledDropDof;
            std::vector<bool> isInteractingBulkCellCenterDof;
            std::vector<bool> isInsideBulkCellCenterDof;
            std::vector<bool> isInteractingBulkFaceDof;
            std::vector<bool> isInsideBulkFaceDof;
            std::vector<bool> isInteractingBulkElement;
            std::vector<bool> isInsideBulkElement;
            std::vector<Scalar> dropOnBulkElementProjectedArea;

            std::unordered_map<std::size_t, std::vector<std::size_t>> dropElementToBulkElementsMap;
            std::unordered_map<std::size_t, std::size_t> bulkElementToDropElementMap;
        } data;



        data.isCoupledDropDof.resize(dropGridGeometry.numDofs(), false);
        data.isInteractingBulkFaceDof.resize(bulkGridGeometry.numFaceDofs(), false);
        data.isInsideBulkFaceDof.resize(bulkGridGeometry.numFaceDofs(), false);
        data.isInteractingBulkCellCenterDof.resize(bulkGridGeometry.numCellCenterDofs(), false);
        data.isInsideBulkCellCenterDof.resize(bulkGridGeometry.numCellCenterDofs(), false);
        data.isInteractingBulkElement.resize(bulkGridGeometry.numCellCenterDofs(), false);
        data.isInsideBulkElement.resize(bulkGridGeometry.numCellCenterDofs(), false);
        data.dropOnBulkElementProjectedArea.resize(bulkGridGeometry.numCellCenterDofs(), 0.0);

        for (const auto& dropElement : elements(dropGridGeometry.gridView()))
        {
             // check if the drop vertices intersect with the bulk grid
            dropFvGeometry.bindElement(dropElement);
            for (const auto& dropScv : scvs(dropFvGeometry))
            {
                //const auto dropVolume = couplingManager.problem(CouplingManager::dropIdx).dropFeature....................................................
                // auto drop = couplingManager.problem(CouplingManager::dropIdx).dropFeatures().droplet(dropScv);
                if (!couplingManager.problem(CouplingManager::dropIdx).dropSolver().hasDroplet(dropElement, dropScv))
                    continue;

                const auto& drop = couplingManager.problem(CouplingManager::dropIdx).dropSolver().droplet(dropElement, dropScv);
                const auto& dropVolume = drop.volume();
                // skip the dof if it is not on the boundary
                // if (dropvolume <1e-15)
                //     continue;

                const auto& dropCenter = drop.center();
                const Scalar dropRadius = drop.radius();

            //    std::cout<<std::endl<<dropRadius<<"   "<<dropCenter<<std::endl;
                const auto& dropElementIdx = drop.elementIdx();
                const auto& dropDofIdx = drop.dofIndex();
                data.isCoupledDropDof[dropDofIdx] = true;

                for (const auto& bulkElement : elements(bulkGridGeometry.gridView()))
                {
                    // check if the drop vertices intersect with the bulk grid
                    bulkFvGeometry.bindElement(bulkElement);
                    const auto bulkElemIdx = bulkGridGeometry.elementMapper().index(bulkElement);

                    const auto coupledFaceDofIndices = coupledFaces(bulkFvGeometry, bulkElement, dropCenter, dropRadius);

                    if(!coupledFaceDofIndices.anyInteractingFace)
                        continue;

                    data.isInteractingBulkCellCenterDof[bulkElemIdx] = true;
                    data.isInteractingBulkElement[bulkElemIdx] = true;
                    bulkGridGeometry.trueDropBoundaryScvf(bulkElemIdx);

                    for (const auto bulkFaceDofIdx : coupledFaceDofIndices.interactingFaces)
                    {
                        data.dropToBulkFaceStencils[dropElementIdx].push_back(bulkFaceDofIdx);
                        data.bulkFaceToDropStencils[bulkFaceDofIdx].push_back(dropDofIdx);
                        data.isInteractingBulkFaceDof[bulkFaceDofIdx] = true;

                    }

                    for (const auto bulkFaceDofIdx : coupledFaceDofIndices.isInsideFace)
                    {
                        data.isInsideBulkFaceDof[bulkFaceDofIdx] = true;
                    }

                    auto elemCenter = bulkElement.geometry().center();
                    auto elemCenterDistance = DropIntersection<Scalar, GlobalPosition>::distancePointToPoint(elemCenter, dropCenter);

                    if (elemCenterDistance < dropRadius)
                    {
                        data.isInsideBulkCellCenterDof[bulkElemIdx] = true;
                   }

                    data.isInsideBulkElement[bulkElemIdx] = coupledFaceDofIndices.allFacesInside;
                    data.dropOnBulkElementProjectedArea[bulkElemIdx] = coupledFaceDofIndices.inFlowProjectedArea;

                    data.dropToBulkCellCenterStencils[dropElementIdx].push_back(bulkElemIdx);
                    data.bulkCellCenterToDropStencils[bulkElemIdx].push_back(dropElementIdx);
                    data.dropElementToBulkElementsMap[dropElementIdx].push_back(bulkElemIdx);
                    data.bulkElementToDropElementMap[bulkElemIdx] = dropElementIdx;

                }
            }
        }

        return data;
    }

//! retruns the relevant coupling maps and stencils
    template<class CouplingManager>
    static void resetStencil(const CouplingManager& couplingManager)
    {

        using std::abs;
        const auto& bulkGridGeometry = couplingManager.problem(CouplingManager::bulkIdx).gridGeometry();
        const auto& dropGridGeometry = couplingManager.problem(CouplingManager::dropIdx).gridGeometry();

        auto dropFvGeometry = localView(dropGridGeometry);
        auto bulkFvGeometry = localView(bulkGridGeometry);

        using BulkGridView = std::decay_t<decltype(bulkGridGeometry.gridView())>;
        using GlobalPosition = typename BulkGridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
        using Scalar = typename BulkGridView::ctype;


                for (const auto& bulkElement : elements(bulkGridGeometry.gridView()))
                {
                    // check if the drop vertices intersect with the bulk grid
                    bulkFvGeometry.bindElement(bulkElement);
                    const auto bulkElemIdx = bulkGridGeometry.elementMapper().index(bulkElement);
                    for (const auto& scvf : scvfs(bulkFvGeometry))
                    {
                        scvf.resetDropBoundary();
                    }
                }
    }

    //! get the bulk faces that are coupled to the drop dofs
    template<class FVElementGeometry, class Element, class GlobalPosition, class Scalar>
    static auto coupledFaces(const FVElementGeometry& fvGeometry,
                              const Element& element,
                              const GlobalPosition& dropCenter,
                              const Scalar& dropRadius)
    {
        struct Result
        {
            std::vector<std::size_t> interactingFaces;
            bool anyInteractingFace;
            bool allFacesInside;
            std::vector<std::size_t> isInsideFace;
            Scalar inFlowProjectedArea;
        } result;

        static constexpr auto dimWorld = GlobalPosition::dimension;

        if(dimWorld == 3)
            coupledFacesThreeDim_(fvGeometry, element, dropCenter, dropRadius, result);
        else
            coupledFacesTwoDim_(fvGeometry, element, dropCenter, dropRadius, result);

        return result;

    }
private:
    //! get the bulk faces that are coupled to the drop dofs
    template<class FVElementGeometry, class Element, class GlobalPosition, class Scalar, class Result>
    static void coupledFacesThreeDim_(const FVElementGeometry& fvGeometry,
                                      const Element& element,
                                      const GlobalPosition& dropCenter,
                                      const Scalar& dropRadius,
                                      Result& result)
    {
        using std::abs;
        static constexpr auto dimWorld = GlobalPosition::dimension;
        auto onDropNormalVector = elemToDropNormalVector_(element, dropCenter);

        // struct Result
        // {
        //     std::vector<std::size_t> interactingFaces;
        //     bool anyInteractingFace;
        //     bool allFacesInside;
        //     std::vector<std::size_t> isInsideFace;
        //     Scalar inFlowProjectedArea;
        // } result;

        result.anyInteractingFace = false;
        result.inFlowProjectedArea = 0.0;
        result.allFacesInside = false;
        int countInsideFaces = 0;
        for (const auto& scvf : scvfs(fvGeometry))
        {
            //scvf.resetDropBoundary();
            auto minDistance = DropIntersection<Scalar, GlobalPosition>::distancePointToSurface(scvf, dropCenter);
            if (minDistance > dropRadius)
                continue;

            result.interactingFaces.push_back(scvf.dofIndex());
            result.anyInteractingFace = true;


            Scalar projectedArea = DropIntersection<Scalar, GlobalPosition>::dropToBulkFaceProjectedArea(dropCenter, dropRadius, scvf);
//std::cout<<std::endl<<scvf.dofIndex()<<"   "<<scvf.area()<<"   "<<projectedArea<<std::endl;
            if (projectedArea > 0.50*scvf.area())
            {
                result.isInsideFace.push_back(scvf.dofIndex());
                //if (scvf.unitOuterNormal()[dimWorld - 1] == 0)
                    scvf.ifDropBoundary();
                    //scvf.setDropInfectedArea(0.0);

                countInsideFaces += 1;
            }
            // else
            // {
            //     scvf.setDropInfectedArea(projectedArea);
            // }


            if (scvf.unitOuterNormal()[0])
            {
                if (onDropNormalVector[0] == 0)
                {
                    Scalar sign = scvf.unitOuterNormal()[0];
                    result.inFlowProjectedArea += sign * projectedArea;
                }
                else
                {
                    Scalar sign = scvf.unitOuterNormal()[0] * onDropNormalVector[0] / std::abs(scvf.unitOuterNormal()[0] * onDropNormalVector[0]);
                    result.inFlowProjectedArea += sign * projectedArea;
                }
            }
        }

        if (countInsideFaces == fvGeometry.numScvf())
            result.allFacesInside = true;

    }

    //! get the bulk faces that are coupled to the drop dofs
    template<class FVElementGeometry, class Element, class GlobalPosition, class Scalar, class Result>
    static void coupledFacesTwoDim_(const FVElementGeometry& fvGeometry,
                              const Element& element,
                              const GlobalPosition& dropCenter,
                              const Scalar& dropRadius,
                              Result& result)
    {
        using std::abs;
        static constexpr auto dimWorld = GlobalPosition::dimension;
        auto onDropNormalVector = elemToDropNormalVector_(element, dropCenter);

        result.anyInteractingFace = false;
        result.inFlowProjectedArea = 0.0;
        result.allFacesInside = false;
        int countInsideFaces = 0;

        for (const auto& scvf : scvfs(fvGeometry))
        {
            auto firstCorner = scvf.corner(0);
            auto secondCorner = scvf.corner(1);
            auto minDistance = DropIntersection<Scalar, GlobalPosition>::distancePointToLine(firstCorner, secondCorner, dropCenter);
            if (minDistance > dropRadius)
                continue;
            result.interactingFaces.push_back(scvf.dofIndex());
            result.anyInteractingFace = true;
            Scalar projectedArea = DropIntersection<Scalar, GlobalPosition>::dropToBulkFaceProjectedArea(dropCenter, dropRadius, scvf);

            // if(scvf.dofIndex() == 35 || scvf.dofIndex() == 42)
            Scalar epsilon = 1e-4;
            if (projectedArea > scvf.area() - std::max(projectedArea, 0.8*scvf.area())*epsilon)
            {
                // std::cout<<"  vvv  "<<scvf.dofIndex()<<std::endl;
                result.isInsideFace.push_back(scvf.dofIndex());
                //if (scvf.unitOuterNormal()[dimWorld - 1] == 0)
                    scvf.setDropBoundary();
                    //scvf.setDropInfectedArea(0.0);

                countInsideFaces += 1;
            }

            if (scvf.unitOuterNormal()[0])
            {
                auto areaExtrusion = DropIntersection<Scalar, GlobalPosition>::dropToBulkFaceProjectedAreaExtrusionFactor(dropCenter, dropRadius, scvf);
                //std::cout<<std::endl<<"HHZZZZZ  "<<scvf.dofIndex()<<"  "<<scvf.area()<<"   "<<areaExtrusion<<std::endl;
                //projectedArea *= areaExtrusion;
                if (onDropNormalVector[0] == 0)
                {
                    Scalar sign = scvf.unitOuterNormal()[0];
                    result.inFlowProjectedArea += sign * areaExtrusion;
                }
                else
                {
                    Scalar sign = scvf.unitOuterNormal()[0] * onDropNormalVector[0] / std::abs(scvf.unitOuterNormal()[0] * onDropNormalVector[0]);
                    result.inFlowProjectedArea += sign * areaExtrusion;
                }
            }
        }
        //std::cout<<std::endl<<"projArea in map  "<<result.inFlowProjectedArea<<std::endl;

        if (countInsideFaces == fvGeometry.numScvf())
            result.allFacesInside = true;
    }

    template<class Element, class GlobalPosition>
    static auto elemToDropNormalVector_(const Element& element,
                                       const GlobalPosition& dropCenter)
    {

        const auto& elementCenter = element.geometry().center();

        auto normalVector = dropCenter;
        normalVector -=  elementCenter;
        normalVector /= normalVector.two_norm();

        return normalVector;
    }

};

} // end namespace Dumux

#endif
