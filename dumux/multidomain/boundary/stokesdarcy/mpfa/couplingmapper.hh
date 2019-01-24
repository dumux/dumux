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
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \copydoc Dumux::StokesDarcyCouplingMapper
 */

#ifndef DUMUX_STOKES_DARCY_COUPLINGMAPPER_HH
#define DUMUX_STOKES_DARCY_COUPLINGMAPPER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <memory>

#include <dune/common/timer.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief Coupling mapper for Stokes and Darcy domains with equal dimension.
 */
template<class MDTraits>
class StokesDarcyCouplingMapper
{
    using Scalar = typename MDTraits::Scalar;

public:
    static constexpr auto stokesCellCenterIdx = typename MDTraits::template DomainIdx<0>();
    static constexpr auto stokesFaceIdx = typename MDTraits::template DomainIdx<1>();
    static constexpr auto cellCenterIdx = typename MDTraits::template DomainIdx<0>();
    static constexpr auto faceIdx = typename MDTraits::template DomainIdx<1>();
    static constexpr auto stokesIdx = stokesCellCenterIdx;
    static constexpr auto darcyIdx = typename MDTraits::template DomainIdx<2>();


private:
    // obtain the type tags of the sub problems
    using StokesTypeTag = typename MDTraits::template SubDomainTypeTag<0>;
    using DarcyTypeTag = typename MDTraits::template SubDomainTypeTag<2>;

    struct ElementMapInfo
    {
        std::size_t eIdx;
        std::size_t scvfIdx;
        std::size_t flipScvfIdx;
    };

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;
    using CouplingManager = typename GET_PROP_TYPE(StokesTypeTag, CouplingManager);

    static_assert(GET_PROP_TYPE(SubDomainTypeTag<stokesIdx>, FVGridGeometry)::discMethod == DiscretizationMethod::staggered,
                  "The free flow domain must use the staggered discretization");

//    static_assert(GET_PROP_TYPE(SubDomainTypeTag<darcyIdx>, FVGridGeometry)::discMethod == DiscretizationMethod::cctpfa,
//                  "The Darcy domain must use the CCTpfa discretization");
public:

    /*!
     * \brief Constructor
     */
    StokesDarcyCouplingMapper(const CouplingManager& couplingManager) : couplingManager_(couplingManager) {}

    /*!
     * \brief Main update routine
     */
//    template<class Stencils>
//    void computeCouplingMapsAndStencils(Stencils& darcyToStokesCellCenterStencils,
//                                        Stencils& darcyToStokesFaceStencils,
//                                        Stencils& stokesCellCenterToDarcyStencils,
//                                        Stencils& stokesFaceToDarcyStencils)
//    {
//        const auto& stokesProblem = couplingManager_.problem(stokesIdx);
//        const auto& darcyProblem = couplingManager_.problem(darcyIdx);
//
//        const auto& stokesFvGridGeometry = stokesProblem.fvGridGeometry();
//        const auto& darcyFvGridGeometry = darcyProblem.fvGridGeometry();
//
//        isCoupledDarcyScvf_.resize(darcyFvGridGeometry.numScvf(), false);
//
//        auto darcyFvGeometry = localView(darcyFvGridGeometry);
//        auto stokesFvGeometry = localView(stokesFvGridGeometry);
//        const auto& stokesGridView = stokesFvGridGeometry.gridView();
//
//        for(const auto& stokesElement : elements(stokesGridView))
//        {
//            stokesFvGeometry.bindElement(stokesElement);
//
//            for(const auto& scvf : scvfs(stokesFvGeometry))
//            {
//                // skip the DOF if it is not on the boundary
//                if(!scvf.boundary())
//                    continue;
//
//                if(scvf.index() == 6)
//                {
//                    Scalar test = 0.0;
//                }
//
//                // get element intersecting with the scvf center
//                // for robustness, add epsilon in unit outer normal direction
//                const auto eps = (scvf.center() - stokesElement.geometry().center()).two_norm()*1e-8;
//                auto globalPos = scvf.center(); globalPos.axpy(eps, scvf.unitOuterNormal());
//                const auto darcyDofIdx = intersectingEntities(globalPos, darcyFvGridGeometry.boundingBoxTree());
//
//                // skip if no intersection was found
//                if(darcyDofIdx.empty())
//                    continue;
//
//                // sanity check
//                if(darcyDofIdx.size() > 1)
//                    DUNE_THROW(Dune::InvalidStateException, "Stokes face dof should only intersect with one Darcy element");
//
//                const auto stokesElementIdx = stokesFvGridGeometry.elementMapper().index(stokesElement);
//                const auto darcyElementIdx = darcyDofIdx[0];
//
//                stokesCellCenterToDarcyStencils[stokesElementIdx].push_back(darcyElementIdx);
//                darcyToStokesCellCenterStencils[darcyElementIdx].push_back(stokesElementIdx);
//                darcyToStokesFaceStencils[darcyElementIdx].push_back(scvf.dofIndex());
//
//                std::cout << "scvf.dofIndex(): " << scvf.dofIndex() << std::endl;
//                std::cout << "scvf.index(): " << scvf.index() << std::endl;
//
//                const auto& darcyElement = darcyFvGridGeometry.element(darcyElementIdx);
//                darcyFvGeometry.bind(darcyElement);
//
//                // find the corresponding Darcy sub control volume face
//                for(const auto& darcyScvf : scvfs(darcyFvGeometry))
//                {
//                    if(!darcyScvf.boundary())
//                        continue;
//
//                    const Scalar distance = (darcyScvf.center() - scvf.center()).two_norm();
//
//                    if(distance < eps)
//                    {
//                        const auto& indexSet = darcyFvGridGeometry.gridInteractionVolumeIndexSets().primaryIndexSet(darcyScvf).nodalIndexSet();
//                        const auto& scvIndices = indexSet.globalScvIndices();
//                        const auto& scvfIndices = indexSet.globalScvfIndices();
//
//                        for(auto&& index : scvIndices)
//                            stokesFaceToDarcyStencils[scvf.dofIndex()].push_back(index);
//
//                        for(auto&& index : scvfIndices)
//                        {
//                            const auto& scvfGlobalDarcy = darcyFvGeometry.scvf(index);
//                            if(!scvfGlobalDarcy.boundary())
//                                continue;
//
//                            for(const auto& data : scvf.pairData())
//                            {
//                                const auto outerParallelFaceDofIdx = data.outerParallelFaceDofIdx;
//                                if(outerParallelFaceDofIdx >= 0)
//                                {
//                                    const auto& scvfGlobalStokes = stokesFvGridGeometry.scvf(data.outerParallelFaceDofIdx);
//                                    if((scvfGlobalDarcy.center() - scvfGlobalStokes.center()).two_norm() < eps)
//                                    {
//                                        darcyToStokesFaceStencils[darcyElementIdx].push_back(scvfGlobalStokes.dofIndex());
//                                    }
//                                }
//
//                            }
//                        }
//
//                        isCoupledDarcyScvf_[darcyScvf.index()] = true;
//                        darcyElementToStokesElementMap_[darcyElementIdx].push_back({stokesElementIdx, scvf.index(), darcyScvf.index()});
//                        stokesElementToDarcyElementMap_[stokesElementIdx].push_back({darcyElementIdx, darcyScvf.index(), scvf.index()});
//                    }
//                }
//            }
//        }
//    }

    template<class Stencils, class StencilsB>
    void computeCouplingMapsAndStencils(Stencils& darcyToStokesCellCenterStencils,
                                        StencilsB& darcyToStokesFaceStencils,
                                        Stencils& stokesCellCenterToDarcyStencils,
                                        Stencils& stokesFaceToDarcyStencils)
    {
        computeCouplingMaps();

        const auto& stokesProblem = couplingManager_.problem(stokesIdx);
        const auto& darcyProblem = couplingManager_.problem(darcyIdx);

        const auto& stokesFvGridGeometry = stokesProblem.fvGridGeometry();
        const auto& darcyFvGridGeometry = darcyProblem.fvGridGeometry();

        auto darcyFvGeometry = localView(darcyFvGridGeometry);
        auto stokesFvGeometry = localView(stokesFvGridGeometry);

        for(const auto& dataHandle : stokesElementToDarcyElementMap_)
        {
            if(dataHandle.second.size() > 1)
                DUNE_THROW(Dune::InvalidStateException, "Stokes face dof should only intersect with one Darcy element");

            const auto& data = dataHandle.second[0];
            const auto stokesElementIdx = dataHandle.first;
            const auto darcyElementIdx = data.eIdx;
            const auto darcyScvfIdx = data.scvfIdx;
            const auto stokesScvfIdx = data.flipScvfIdx;
            const auto& stokesScvf = stokesFvGridGeometry.scvf(stokesScvfIdx);

            const auto& darcyElement = darcyFvGridGeometry.element(darcyElementIdx);
            darcyFvGeometry.bind(darcyElement);
            const auto& darcyScvf = darcyFvGeometry.scvf(darcyScvfIdx);

            stokesCellCenterToDarcyStencils[stokesElementIdx].push_back(darcyElementIdx);
            darcyToStokesCellCenterStencils[darcyElementIdx].push_back(stokesElementIdx);
            darcyToStokesFaceStencils[darcyElementIdx].first.push_back(stokesScvf.dofIndex());
            darcyToStokesFaceStencils[darcyElementIdx].second.push_back(stokesScvf.index());

            const auto& indexSet = darcyFvGridGeometry.gridInteractionVolumeIndexSets().primaryIndexSet(darcyScvf).nodalIndexSet();
            const auto& scvIndices = indexSet.globalScvIndices();
            const auto& scvfIndices = indexSet.globalScvfIndices();

            for(auto&& index : scvIndices)
            {
                stokesFaceToDarcyStencils[stokesScvf.dofIndex()].push_back(index);
                darcyToStokesFaceStencils[index].first.push_back(stokesScvf.dofIndex());
                darcyToStokesFaceStencils[index].second.push_back(stokesScvf.index());
                darcyElementToStokesElementMap_[index].push_back({stokesElementIdx, stokesScvfIdx, darcyScvfIdx});
            }

            for(auto&& index : scvfIndices)
            {
                const auto& scvfIVDarcy = darcyFvGeometry.scvf(index);
                if(!scvfIVDarcy.boundary())
                    continue;

                if(!darcyElementToStokesElementMap_.count(scvfIVDarcy.insideScvIdx()))
                    continue;

                // prepare the coupling context
                const auto& stokesElementIndices = darcyElementToStokesElementMap_.at(scvfIVDarcy.insideScvIdx());

                for(auto&& indices : stokesElementIndices)
                {
                    if(indices.flipScvfIdx == index)
                    {
                        darcyToStokesFaceStencils[darcyElementIdx].first.push_back(stokesFvGridGeometry.scvf(indices.scvfIdx).dofIndex());
                        darcyToStokesFaceStencils[darcyElementIdx].second.push_back(stokesFvGridGeometry.scvf(indices.scvfIdx).index());
                    }
                }
            }

        }
    }

    void computeCouplingMaps()
    {
        const auto& stokesProblem = couplingManager_.problem(stokesIdx);
        const auto& darcyProblem = couplingManager_.problem(darcyIdx);

        const auto& stokesFvGridGeometry = stokesProblem.fvGridGeometry();
        const auto& darcyFvGridGeometry = darcyProblem.fvGridGeometry();

        isCoupledDarcyScvf_.resize(darcyFvGridGeometry.numScvf(), false);

        auto darcyFvGeometry = localView(darcyFvGridGeometry);
        auto stokesFvGeometry = localView(stokesFvGridGeometry);
        const auto& stokesGridView = stokesFvGridGeometry.gridView();

        for(const auto& stokesElement : elements(stokesGridView))
        {
            stokesFvGeometry.bindElement(stokesElement);

            for(const auto& scvf : scvfs(stokesFvGeometry))
            {
                // skip the DOF if it is not on the boundary
                if(!scvf.boundary())
                    continue;

                // get element intersecting with the scvf center
                // for robustness, add epsilon in unit outer normal direction
                const auto eps = (scvf.center() - stokesElement.geometry().center()).two_norm()*1e-8;
                auto globalPos = scvf.center(); globalPos.axpy(eps, scvf.unitOuterNormal());
                const auto darcyDofIdx = intersectingEntities(globalPos, darcyFvGridGeometry.boundingBoxTree());

                // skip if no intersection was found
                if(darcyDofIdx.empty())
                    continue;

                // sanity check
                if(darcyDofIdx.size() > 1)
                    DUNE_THROW(Dune::InvalidStateException, "Stokes face dof should only intersect with one Darcy element");

                const auto stokesElementIdx = stokesFvGridGeometry.elementMapper().index(stokesElement);
                const auto darcyElementIdx = darcyDofIdx[0];

                const auto& darcyElement = darcyFvGridGeometry.element(darcyElementIdx);
                darcyFvGeometry.bindElement(darcyElement);

                // find the corresponding Darcy sub control volume face
                for(const auto& darcyScvf : scvfs(darcyFvGeometry))
                {
                    if(!darcyScvf.boundary())
                        continue;

                    const Scalar distance = (darcyScvf.center() - scvf.center()).two_norm();

                    if(distance < eps)
                    {
                        isCoupledDarcyScvf_[darcyScvf.index()] = true;
                        darcyElementToStokesElementMap_[darcyElementIdx].push_back({stokesElementIdx, scvf.index(), darcyScvf.index()});
                        stokesElementToDarcyElementMap_[stokesElementIdx].push_back({darcyElementIdx, darcyScvf.index(), scvf.index()});
                    }
                }
            }
        }
    }

    /*!
     * \brief Returns whether a Darcy scvf is coupled to the other domain
     */
    bool isCoupledDarcyScvf(std::size_t darcyScvfIdx) const
    {
        return isCoupledDarcyScvf_[darcyScvfIdx];
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
    std::unordered_map<std::size_t, std::vector<ElementMapInfo>> darcyElementToStokesElementMap_;
    std::unordered_map<std::size_t, std::vector<ElementMapInfo>> stokesElementToDarcyElementMap_;

    std::vector<bool> isCoupledDarcyScvf_;

    const CouplingManager& couplingManager_;
};

} // end namespace Dumux

#endif
