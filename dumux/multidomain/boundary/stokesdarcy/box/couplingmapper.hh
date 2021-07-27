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
#include <dumux/geometry/geometryintersection.hh>
#include <dumux/geometry/intersectingentities.hh>

#include <dumux/discretization/method.hh>
#include <dumux/freeflow/slipcondition.hh>

namespace Dumux {

template<class CouplingFacet, class Container, class IndexSet, class Iterator>
class FacetIterator : public Dune::ForwardIteratorFacade<FacetIterator<CouplingFacet,
                                                                       Container,
                                                                       IndexSet,
                                                                       Iterator>,
                                                         const CouplingFacet>
{
    using ThisType = FacetIterator<CouplingFacet, Container, IndexSet, Iterator>;
public:
    FacetIterator(const Iterator& it, const Container& container, const IndexSet& indices)
    : it_(it), container_(&container), indices_(&indices) {}

    FacetIterator() : it_(Iterator()), container_(nullptr), indices_(nullptr) {}

    //! dereferencing yields a coupling facet
    const CouplingFacet& dereference() const
    {
        return container_->at(*it_);
    }

    bool equals(const ThisType& other) const
    {
        return it_ == other.it_;
    }

    template<bool enable = std::is_integral<IndexSet>::value, std::enable_if_t<!enable, int> = 0>
    void increment()
    {
        it_++;
        if(it_ == (*indices_)[localIdx_].end())
        {
            localIdx_++;
            auto it = std::find_if (indices_->begin()+localIdx_, indices_->end(), [] (const auto& idx) { return idx.size() > 0; });

            if(it != indices_->end())
                it_ = it->begin();
        }
    }

    template<bool enable = std::is_integral<IndexSet>::value, std::enable_if_t<enable, int> = 0>
    void increment()
    {
        it_++;
    }

private:
    Iterator it_;
    const Container* container_;
    const IndexSet* indices_;
    std::size_t localIdx_{0};
};

template<class Mapper, std::size_t i>
auto couplingFacets(Dune::index_constant<i> domainI, const Mapper& mapper, std::size_t eIdx)
{
    using FacetContainer = typename Mapper::FacetContainer;
    using IndexContainerScvf = typename Mapper::IndexContainerScvf;
    using IndexContainerElement = typename Mapper::IndexContainerElement;
    using FacetIterator = FacetIterator<typename FacetContainer::value_type, FacetContainer, IndexContainerElement, typename IndexContainerScvf::const_iterator>;
    const auto& indexSet = mapper.couplingFacetIdxMap(domainI).at(eIdx);
    auto itBegin = std::find_if (indexSet.begin(), indexSet.end(), [] (const auto& idx) { return idx.size() > 0; });
    auto itEnd = std::find_if(std::reverse_iterator(indexSet.end()),
                              std::reverse_iterator(indexSet.begin()),
                                [] (const auto& idx) { return idx.size() > 0; });

    return Dune::IteratorRange<FacetIterator>(FacetIterator(itBegin->begin(), mapper.couplingFacets(), indexSet),
                                              FacetIterator(itEnd->end(), mapper.couplingFacets(), indexSet));
}

template<class Mapper, std::size_t i>
auto couplingFacets(Dune::index_constant<i> domainI, const Mapper& mapper, std::size_t eIdx, std::size_t localScvfIdx)
{
    using FacetContainer = typename Mapper::FacetContainer;
    using IndexContainerScvf = typename Mapper::IndexContainerScvf;
    using FacetIterator = FacetIterator<typename FacetContainer::value_type, FacetContainer, IndexContainerScvf, typename IndexContainerScvf::const_iterator>;
    const auto& indexSet = mapper.couplingFacetIdxMap(domainI).at(eIdx)[localScvfIdx];
    return Dune::IteratorRange<FacetIterator>(FacetIterator(indexSet->begin(), mapper.couplingFacets(), indexSet),
                                              FacetIterator(indexSet->end(), mapper.couplingFacets(), indexSet));
}

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling mapper for Stokes and Darcy domains with equal dimension when using the Box scheme for the Darcy domain.
 */
template<class MDTraits>
class StokesDarcyCouplingMapperBox
{
private:
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using Scvf = typename GridGeometry<id>::SubControlVolumeFace::Traits::Geometry;

    static constexpr auto freeFlowIdx = FreeFlowPorousMediumCouplingManagerBase<MDTraits>::freeFlowIdx;
    static constexpr auto porousMediumIdx = FreeFlowPorousMediumCouplingManagerBase<MDTraits>::porousMediumIdx;

    static constexpr int dim = GridView<porousMediumIdx>::dimension;
    static constexpr int dimWorld = GridView<porousMediumIdx>::dimensionworld;

    using ctype = typename Dune::PromotionTraits<typename GridView<freeFlowIdx>::ctype,
                                                 typename GridView<porousMediumIdx>::ctype>::PromotedType;

    static_assert(GridView<freeFlowIdx>::dimension == dim, "The grids must have the same dimension");
    static_assert(GridView<freeFlowIdx>::dimensionworld == dimWorld, "The grids must have the same world dimension");
    static_assert(GridGeometry<freeFlowIdx>::discMethod == DiscretizationMethod::staggered, "The free flow domain must use the staggered discretization");
    static_assert(GridGeometry<porousMediumIdx>::discMethod == DiscretizationMethod::box, "The Darcy domain must use the Box discretization");

public:
    // export the type describing a coupling segment
    struct CouplingFacet
    {
        // each intersection segment is described by a simplex geometry of codimension one
        using Geometry = Dune::AffineGeometry<ctype, dim-1, dimWorld>;

        std::size_t ffEIdx;
        std::size_t pmEIdx;
        std::size_t ffScvfIdx;
        std::size_t pmScvfIdx;
        std::size_t idx;
        Geometry geometry;
    };

    using FacetContainer = std::vector<CouplingFacet>;
    using IndexContainerScvf = std::vector<std::size_t>;
    using IndexContainerElement = std::vector<IndexContainerScvf>;

    // template<class IndexSet, bool enable = std::is_integral<IndexSet>::value, std::enable_if_t<enable, int> = 0>
    // auto couplingFacets(const std::vector<IndexSet>& indexSet)
    // {
    //     using FacetIterator = FacetIterator<CouplingFacet, FacetContainer, std::vector<IndexSet>, typename std::vector<IndexSet>::const_iterator>;
    //     return Dune::IteratorRange<FacetIterator>(FacetIterator(indexSet->begin(), couplingFacets_, indexSet),
    //                                               FacetIterator(indexSet->end(), couplingFacets_, indexSet));
    // }

    // template<class IndexSet, bool enable = std::is_integral<IndexSet>::value, std::enable_if_t<!enable, int> = 0>
    // auto couplingFacets(const std::vector<IndexSet>& indexSet)
    // {
    //     using FacetIterator = FacetIterator<CouplingFacet, FacetContainer, std::vector<IndexSet>, typename IndexSet::const_iterator>;
    //     auto itBegin = std::find_if (indexSet.begin(), indexSet.end(), [] (const auto& idx) { return idx.size() > 0; });
    //     auto itEnd = std::find_if(std::reverse_iterator(indexSet.end()),
    //                               std::reverse_iterator(indexSet.begin()),
    //                               [] (const auto& idx) { return idx.size() > 0; });

    //     return Dune::IteratorRange<FacetIterator>(FacetIterator(itBegin->begin(), couplingFacets_, indexSet),
    //                                               FacetIterator(itEnd->end(), couplingFacets_, indexSet));
    // }

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

        for(const auto& elementData : ffCouplingFacetIdxMap_)
        {
            const auto stokesEIdx = elementData.first;

            for(const auto& couplingFacet : Dumux::couplingFacets(Dune::index_constant<freeFlowIdx>(), *this, stokesEIdx))
            {
                const auto darcyEIdx = couplingFacet.pmEIdx;
                const auto stokesScvfIdx = couplingFacet.ffScvfIdx;
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

                if(!(slipCondition() == SlipCondition::BJS))
                {
                    const std::size_t numSubFaces = stokesScvf.pairData().size();
                    // Account for all interior sub-faces which include data from a boundary with slip condition
                    for (int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
                    {
                        const auto eIdx = stokesScvf.insideScvIdx();
                        const auto& lateralStokesScvf = stokesFvGeometry.scvf(eIdx, stokesScvf.pairData(localSubFaceIdx).localLateralFaceIdx);
                        if(lateralStokesScvf.dofIndex() != stokesScvf.dofIndex() && !lateralStokesScvf.boundary())
                            for (auto&& scv : scvs(darcyFvGeometry))
                                stokesFaceToDarcyStencils[lateralStokesScvf.dofIndex()].push_back(scv.dofIndex());
                    }
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

        std::size_t couplingFaceIdx = 0;
        // find all darcy faces coupling to stokes
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
                        using IntersectionAlgorithm = GeometryIntersection<Scvf<porousMediumIdx>, Scvf<freeFlowIdx>>;
                        typename IntersectionAlgorithm::Intersection rawIs;
                        if(IntersectionAlgorithm::intersection(darcyScvfGeometry, stokesScvf.geometry(), rawIs))
                        {
                            if(pmCouplingFacetIdxMap_[darcyEIdx].size() == 0)
                                pmCouplingFacetIdxMap_[darcyEIdx].resize(darcyFvGeometry.numScvf());

                            if(ffCouplingFacetIdxMap_[stokesEIdx].size() == 0)
                                ffCouplingFacetIdxMap_[stokesEIdx].resize(stokesFvGeometry.numScvf());

                            const auto is = typename CouplingFacet::Geometry(Dune::GeometryTypes::simplex(dim-1), rawIs);
                            pmCouplingFacetIdxMap_[darcyEIdx][darcyScvf.index()].push_back(couplingFaceIdx);
                            ffCouplingFacetIdxMap_[stokesEIdx][stokesScvf.localFaceIdx()].push_back(couplingFaceIdx);
                            couplingFacets_.push_back({stokesEIdx, darcyEIdx, stokesScvf.index(), darcyScvf.index(), couplingFaceIdx, is});
                            couplingFaceIdx++;
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
        return (pmCouplingFacetIdxMap_.count(eIdx) == 0) ?  false
                                                         : pmCouplingFacetIdxMap_.at(eIdx)[scvfLocalIdx].size() > 0;
    }

    /*!
     * \brief A map that returns all elements coupled to the other domain
     */
    template<std::size_t i>
    const auto& couplingFacetIdxMap(Dune::index_constant<i> domainI) const
    {
        if constexpr (i == porousMediumIdx)
            return pmCouplingFacetIdxMap_;
        else if constexpr (i == freeFlowIdx)
            return ffCouplingFacetIdxMap_;
    }

    const CouplingFacet& couplingFacet(std::size_t idx) const
    {
        return couplingFacets_[idx];
    }

    const FacetContainer& couplingFacets() const
    {
        return couplingFacets_;
    }


private:
    FacetContainer couplingFacets_;
    std::unordered_map<std::size_t, IndexContainerElement> pmCouplingFacetIdxMap_;
    std::unordered_map<std::size_t, IndexContainerElement> ffCouplingFacetIdxMap_;
};

} // end namespace Dumux

#endif
