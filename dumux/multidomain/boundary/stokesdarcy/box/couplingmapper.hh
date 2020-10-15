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
 * \ingroup StokesDarcyCoupling
 * \copydoc Dumux::StokesDarcyCouplingMapper
 */

#ifndef DUMUX_STOKES_DARCY_COUPLINGMAPPER_BOX_HH
#define DUMUX_STOKES_DARCY_COUPLINGMAPPER_BOX_HH

#include <type_traits>
#include <unordered_map>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/promotiontraits.hh>
#include <dune/geometry/affinegeometry.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/geometry/geometryintersection.hh>
#include <dumux/common/geometry/intersectingentities.hh>

#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling mapper for Stokes and Darcy domains with equal dimension when using the Box scheme for the Darcy domain.
 */
template<class MDTraits>
class StokesDarcyCouplingMapperBox
{
private:
    // obtain the type tags of the sub problems
    using StokesTypeTag = typename MDTraits::template SubDomain<0>::TypeTag;
    using DarcyTypeTag = typename MDTraits::template SubDomain<2>::TypeTag;

    // sub domain grid geometries & scvf geometries
    using StokesGG = GetPropType<StokesTypeTag, Properties::GridGeometry>;
    using DarcyGG = GetPropType<DarcyTypeTag, Properties::GridGeometry>;
    using StokesScvfGeometry = typename StokesGG::SubControlVolumeFace::Traits::Geometry;
    using DarcyScvfGeometry = typename DarcyGG::SubControlVolumeFace::Traits::Geometry;
    using ctype = typename Dune::PromotionTraits<typename StokesGG::GridView::ctype,
                                                 typename DarcyGG::GridView::ctype>::PromotedType;

    static constexpr int dim = DarcyGG::GridView::dimension;
    static constexpr int dimWorld = DarcyGG::GridView::dimensionworld;

    static_assert(StokesGG::GridView::dimension == dim, "The grids must have the same dimension");
    static_assert(StokesGG::GridView::dimensionworld == dimWorld, "The grids must have the same world dimension");
    static_assert(StokesGG::discMethod == DiscretizationMethod::staggered, "The free flow domain must use the staggered discretization");
    static_assert(DarcyGG::discMethod == DiscretizationMethod::box, "The Darcy domain must use the Box discretization");

public:

    // export the type describing a coupling segment
    struct CouplingSegment
    {
        // each intersection segment is described by a simplex geometry of codimension one
        using Geometry = Dune::AffineGeometry<ctype, dim-1, dimWorld>;

        std::size_t eIdx;
        std::size_t scvfIdx;
        std::size_t flipScvfIdx;
        Geometry geometry;
    };

    /*!
     * \brief Main update routine
     */
    template<class CouplingManager, class StencilA, class StencilB>
    void computeCouplingMapsAndStencils(const CouplingManager& couplingManager,
                                        StencilA& darcyToStokesCellCenterStencils,
                                        StencilB& darcyToStokesFaceStencils,
                                        StencilA& stokesCellCenterToDarcyStencils,
                                        StencilA& stokesFaceToDarcyStencils)
    {
        computeCouplingMaps(couplingManager);

        const auto& stokesProblem = couplingManager.problem(CouplingManager::freeFlowIdx);
        const auto& darcyProblem = couplingManager.problem(CouplingManager::porousMediumIdx);

        const auto& stokesFvGridGeometry = stokesProblem.gridGeometry();
        const auto& darcyFvGridGeometry = darcyProblem.gridGeometry();

        for(const auto& dataHandle : stokesElementToDarcyElementMap_)
        {
            const auto stokesEIdx = dataHandle.first;

            for (const auto& darcyData : dataHandle.second)
            {
                const auto darcyEIdx = darcyData.eIdx;
                const auto stokesScvfIdx = darcyData.flipScvfIdx;
                const auto& stokesScvf = stokesFvGridGeometry.scvf(stokesScvfIdx);

                const auto& darcyElement = darcyFvGridGeometry.element(darcyEIdx);
                auto darcyFvGeometry = localView(darcyFvGridGeometry);
                darcyFvGeometry.bind(darcyElement);

                const auto stokesElement = stokesFvGridGeometry.element(stokesEIdx);
                auto stokesFvGeometry = localView(stokesFvGridGeometry);
                stokesFvGeometry.bind(stokesElement);

                darcyToStokesCellCenterStencils[darcyEIdx].push_back(stokesEIdx);
                darcyToStokesFaceStencils[darcyEIdx].first.push_back(stokesScvf.dofIndex());
                darcyToStokesFaceStencils[darcyEIdx].second.push_back(stokesScvf.index());

                for (auto&& scv : scvs(darcyFvGeometry))
                {
                    stokesCellCenterToDarcyStencils[stokesEIdx].push_back(scv.dofIndex());
                    stokesFaceToDarcyStencils[stokesScvf.dofIndex()].push_back(scv.dofIndex());
                }

                const std::size_t numSubFaces = stokesScvf.pairData().size();

                // Account for all sub faces. This is needed when a slip condition is set.
                for (int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
                {
                    const auto eIdx = stokesScvf.insideScvIdx();
                    // Get the face normal to the face the dof lives on. The staggered sub face conincides with half of this lateral face.
                    const auto& lateralStokesScvf = stokesFvGeometry.scvf(eIdx, stokesScvf.pairData(localSubFaceIdx).localLateralFaceIdx);
                    for (auto&& scv : scvs(darcyFvGeometry))
                        if(lateralStokesScvf.dofIndex() != stokesScvf.dofIndex() && !lateralStokesScvf.boundary())
                            stokesFaceToDarcyStencils[lateralStokesScvf.dofIndex()].push_back(scv.dofIndex());
                }

            }
        }
    }

    template<class CouplingManager>
    void computeCouplingMaps(const CouplingManager& couplingManager)
    {
        const auto& stokesProblem = couplingManager.problem(CouplingManager::freeFlowIdx);
        const auto& darcyProblem = couplingManager.problem(CouplingManager::porousMediumIdx);

        const auto& stokesFvGridGeometry = stokesProblem.gridGeometry();
        const auto& darcyFvGridGeometry = darcyProblem.gridGeometry();

        // find all darcy faces coupling to stokes
        isCoupledDarcyScvf_.resize(darcyFvGridGeometry.gridView().size(0));
        for (const auto& darcyElement : elements(darcyFvGridGeometry.gridView()))
        {
            const auto darcyEIdx = darcyFvGridGeometry.elementMapper().index(darcyElement);
            auto darcyFvGeometry = localView(darcyFvGridGeometry);
            darcyFvGeometry.bindElement(darcyElement);

            for (const auto& darcyScvf : scvfs(darcyFvGeometry))
            {
                if (!darcyScvf.boundary())
                    continue;

                // find all stokes elements that intersect with the face
                const auto& darcyScvfGeometry = darcyScvf.geometry();
                const auto rawIntersections = intersectingEntities(darcyScvfGeometry, stokesFvGridGeometry.boundingBoxTree());
                if (rawIntersections.empty())
                    continue;

                if (isCoupledDarcyScvf_[darcyEIdx].empty())
                    isCoupledDarcyScvf_[darcyEIdx].assign(darcyFvGeometry.numScvf(), false);

                for (const auto& rawIntersection : rawIntersections)
                {
                    const auto stokesEIdx = rawIntersection.second();
                    const auto stokesElement = stokesFvGridGeometry.element(stokesEIdx);
                    auto stokesFvGeometry = localView(stokesFvGridGeometry);
                    stokesFvGeometry.bindElement(stokesElement);

                    for (const auto& stokesScvf : scvfs(stokesFvGeometry))
                    {
                        if (!stokesScvf.boundary())
                            continue;

                        // intersect the geometries
                        using IntersectionAlgorithm = GeometryIntersection<DarcyScvfGeometry, StokesScvfGeometry>;
                        typename IntersectionAlgorithm::Intersection rawIs;
                        if(IntersectionAlgorithm::intersection(darcyScvfGeometry, stokesScvf.geometry(), rawIs))
                        {
                            const auto is = typename CouplingSegment::Geometry(Dune::GeometryTypes::simplex(dim-1), rawIs);
                            isCoupledDarcyScvf_[darcyEIdx][darcyScvf.index()] = true;
                            darcyElementToStokesElementMap_[darcyEIdx].push_back({stokesEIdx, stokesScvf.index(), darcyScvf.index(), is});
                            stokesElementToDarcyElementMap_[stokesEIdx].push_back({darcyEIdx, darcyScvf.index(), stokesScvf.index(), is});
                        }
                    }
                }
            }
        }
    }

    /*!
     * \brief Returns whether a Darcy scvf is coupled to the other domain
     */
    bool isCoupledDarcyScvf(std::size_t eIdx, std::size_t scvfLocalIdx) const
    {
        if(isCoupledDarcyScvf_[eIdx].size() > 0)
            return isCoupledDarcyScvf_[eIdx][scvfLocalIdx];

        return false;
    }


    /*!
     * \brief A map that returns all Stokes elements coupled to a Darcy element
     */
    const auto& darcyElementToStokesElementMap() const
    {
        return darcyElementToStokesElementMap_;
    }

    /*!
     * \brief A map that returns all Darcy elements coupled to a Stokes element
     */
    const auto& stokesElementToDarcyElementMap() const
    {
        return stokesElementToDarcyElementMap_;
    }

private:
    std::unordered_map<std::size_t, std::vector<CouplingSegment>> darcyElementToStokesElementMap_;
    std::unordered_map<std::size_t, std::vector<CouplingSegment>> stokesElementToDarcyElementMap_;

    std::vector<std::vector<bool>> isCoupledDarcyScvf_;
};

} // end namespace Dumux

#endif
