// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup DualNetworkCoupling
 * \ingroup PoreNetworkModels
 * \copydoc Dumux::PoreNetwork::DualNetworkCouplingMapper
 */

#ifndef DUMUX_DUAL_NETWORK_COUPLINGMAPPER_HH
#define DUMUX_DUAL_NETWORK_COUPLINGMAPPER_HH

#include <type_traits>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>

#include <dune/common/iteratorrange.hh>
#include <dune/common/iteratorfacades.hh>

#include <dune/common/exceptions.hh>
#include <dumux/common/entitymap.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup DualNetworkCoupling
 * \ingroup PoreNetworkModels
 * \brief Coupling mapper for Stokes and Darcy domains with equal dimension.
 */

template<class Scalar>
class DualNetworkCouplingMapper
{
    struct HostGridConnectionInfo
    {
        std::size_t hostGridElementIndex;
        std::size_t voidVertexHostIdx;
        std::size_t solidVertexHostIdx;
        Scalar connectionArea;
        Scalar connectionLength;
        std::vector<std::size_t> voidElementHostIdx;
        std::vector<std::size_t> solidElementHostIdx;
        std::vector<std::size_t> coupledVoidVertexHostIdx;
        std::vector<std::size_t> coupledSolidVertexHostIdx;
        std::size_t connectionGlobalId;
    };

    struct SubGridConnectionInfo
    {
        std::size_t id; // the global id of the connections
        std::size_t solidVertexIdx; // the directly coupled solid vertex
        std::size_t voidVertexIdx; // the directly coupled void vertex
        std::size_t someSolidElementIdx; // index of one the solid elements adjacent to the solid vertex
        std::size_t someVoidElementIdx; // index of one the void elements adjacent to the solid vertex
        std::vector<std::size_t> convectionVoidElementIdx; // all void elements adjacent to the own void vertex coupled to the same solid vertex
        Scalar connectionArea;
        Scalar connectionLength;
    };

    template<class Vector>
    class ConnectionIterator : public Dune::ForwardIteratorFacade<ConnectionIterator<Vector>, const SubGridConnectionInfo>
    {
        using ThisType = ConnectionIterator<Vector>;
        using Iterator = typename Vector::const_iterator;
    public:
        ConnectionIterator(const Iterator& it, const std::vector<SubGridConnectionInfo>& info)
        : it_(it), InfoPtr_(&info) {}

        ConnectionIterator() : it_(Iterator()), InfoPtr_(nullptr) {}

        //! dereferencing yields a subcontrol volume
        const SubGridConnectionInfo& dereference() const
        { return InfoPtr_->at(*it_); }
        // { return (*InfoPtr_)[*it_]; }

        bool equals(const ThisType& other) const
        { return it_ == other.it_; }

        void increment()
        { ++it_; }

    private:
        Iterator it_;
        const std::vector<SubGridConnectionInfo>* InfoPtr_;
    };

public:
    using Stencil = std::vector<std::size_t>;

    template<class HostGridView, class HostGridData, class VoidGridGeometry, class SolidGridGeometry>
    DualNetworkCouplingMapper(const HostGridView& hostGridView,
                              const HostGridData& hostGridData,
                              const VoidGridGeometry& voidGridGeometry,
                              const SolidGridGeometry& solidGridGeometry)
    {
        fillVertexMap_(hostGridView, voidGridGeometry, voidHostToSubVertexIdxMap_);
        fillVertexMap_(hostGridView, solidGridGeometry, solidHostToSubVertexIdxMap_);
        fillElementMap_(hostGridView, voidGridGeometry, voidHostToSubElementIdxMap_);
        fillElementMap_(hostGridView, solidGridGeometry, solidHostToSubElementIdxMap_);

        isCoupledVoidDof_.resize(voidGridGeometry.gridView().size(1), false);
        isCoupledSolidDof_.resize(solidGridGeometry.gridView().size(1), false);

        const auto connectionInfo = getConnectionInfo_(hostGridView, hostGridData);
        connectionInfo_.resize(connectionInfo.size());
        for (const auto& info : connectionInfo)
        {
            auto voidHostToSubVertexIdx = [&](const auto hostIdx)
            { return voidHostToSubVertexIdxMap_.at(hostIdx); };

            auto solidHostToSubVertexIdx = [&](const auto hostIdx)
            { return solidHostToSubVertexIdxMap_.at(hostIdx); };

            const auto directlyCoupledVoidDofIdx = voidHostToSubVertexIdx(info.voidVertexHostIdx);
            const auto directlyCoupledSolidDofIdx = solidHostToSubVertexIdx(info.solidVertexHostIdx);

            auto coupledVoidElementIdxSub = info.voidElementHostIdx;
            auto coupledSolidElementIdxSub = info.solidElementHostIdx;

            // convert hostgrid indices to subgrid indices
            std::transform(coupledVoidElementIdxSub.begin(), coupledVoidElementIdxSub.end(),
                           coupledVoidElementIdxSub.begin(), [&](const auto eIdx){ return voidHostToSubElementIdxMap_.at(eIdx); });
            std::transform(coupledSolidElementIdxSub.begin(), coupledSolidElementIdxSub.end(),
                           coupledSolidElementIdxSub.begin(), [&](const auto eIdx){ return solidHostToSubElementIdxMap_.at(eIdx); });

            // initialize an empty vector - will be filled later in a second loop
            auto convectionVoidElementIdx = std::vector<std::size_t>();

            voidToSolidConnectionIds_[directlyCoupledVoidDofIdx].emplace_back(info.connectionGlobalId);
            solidToVoidConnectionIds_[directlyCoupledSolidDofIdx].emplace_back(info.connectionGlobalId);

            connectionInfo_[info.connectionGlobalId] = SubGridConnectionInfo{info.connectionGlobalId,
                                                                             directlyCoupledSolidDofIdx,
                                                                             directlyCoupledVoidDofIdx,
                                                                             coupledSolidElementIdxSub[0],
                                                                             coupledVoidElementIdxSub[0],
                                                                             convectionVoidElementIdx,
                                                                             info.connectionArea,
                                                                             info.connectionLength};

            hostGridElementIndexToGlobalId_[info.hostGridElementIndex] = info.connectionGlobalId;

            isCoupledVoidDof_[directlyCoupledVoidDofIdx] = true;
            isCoupledSolidDof_[directlyCoupledSolidDofIdx] = true;

            for (const auto eIdxVoidHost : info.voidElementHostIdx)
            {
                const auto eIdxSubVoid = voidHostToSubElementIdxMap_.at(eIdxVoidHost);
                voidToSolidStencils_[eIdxSubVoid].push_back(directlyCoupledSolidDofIdx);
            }

            for (const auto eIdxSolidHost : info.solidElementHostIdx)
            {
                const auto eIdxSubSolid = solidHostToSubElementIdxMap_.at(eIdxSolidHost);
                solidToVoidStencils_[eIdxSubSolid].push_back(directlyCoupledVoidDofIdx);
            }
        }

        for (auto& stencil : voidToSolidStencils_)
            removeDuplicates_(stencil.second);
        for (auto& stencil : solidToVoidStencils_)
            removeDuplicates_(stencil.second);

        // second loop: find void element coupling for convective transport
        auto voidFVGeometry = localView(voidGridGeometry);
        for (const auto& voidElement : elements(voidGridGeometry.gridView()))
        {
            voidFVGeometry.bindElement(voidElement);
            std::array<std::size_t, 2> dofIndex;
            std::array<std::vector<std::size_t>, 2> coupledSolidVertexIdx;
            for (const auto& scv : scvs(voidFVGeometry))
                dofIndex[scv.indexInElement()] = scv.dofIndex();

            // check if both void pores couple to the same solid pore
            if (isCoupledVoidDof_[dofIndex[0]] && isCoupledVoidDof_[dofIndex[1]])
            {
                for (auto& conn0 : voidToSolidConnectionIds_[dofIndex[0]])
                {
                    for (auto& conn1 : voidToSolidConnectionIds_[dofIndex[1]])
                    {
                        const auto globalId0 = conn0;
                        const auto globalId1 = conn1;
                        assert(globalId0 != globalId1);

                        if (connectionInfo_[globalId0].solidVertexIdx == connectionInfo_[globalId1].solidVertexIdx)
                        {
                            const auto voidElemIdx = voidGridGeometry.elementMapper().index(voidElement);
                            connectionInfo_[globalId0].convectionVoidElementIdx.push_back(voidElemIdx);
                            connectionInfo_[globalId1].convectionVoidElementIdx.push_back(voidElemIdx);
                        }
                    }
                }
            }
        }

        for (auto& entry : voidToSolidConnectionIds_)
        {
            removeDuplicates_(entry.second);

            std::cout << "void dof " << entry.first << " couples to " << entry.second.size() << " solid dofs: " << std::endl;
            for (auto& conn : entry.second)
            {
                const auto& info = connectionInfo_[conn];
                assert(entry.first == info.voidVertexIdx);
                std::cout << "solid vertex " << info.solidVertexIdx << " with elems ";
                for (const auto e : info.convectionVoidElementIdx)
                    std::cout << e << " ";
                std:: cout << "||" << std::endl;
            }

            std::cout << std::endl;
        }

        for (auto& entry : solidToVoidConnectionIds_)
        {
            removeDuplicates_(entry.second);

            std::cout << "solid dof " << entry.first << " couples to " << entry.second.size() << " void dofs: " << std::endl;
            for (auto& conn : entry.second)
            {
                const auto& info = connectionInfo_[conn];
                assert(entry.first == info.solidVertexIdx);
                std::cout << "void vertex " << info.voidVertexIdx << " with elems ";
                for (const auto e : info.convectionVoidElementIdx)
                    std::cout << e << " ";
                std:: cout << "||" << std::endl;
            }

            std::cout << std::endl;
        }

        // TODO maybe delete hostToSub maps? or make public?
    }

    const auto& voidToSolidStencils() const
    { return voidToSolidStencils_; }

    const auto& solidToVoidStencils() const
    { return solidToVoidStencils_; }

    const std::vector<bool>& isCoupledVoidDof() const
    { return isCoupledVoidDof_; }

    const std::vector<bool>& isCoupledSolidDof() const
    { return isCoupledSolidDof_; }

    //! Returns an iterator allowing for (const auto& conn : voidToSolidConnections(dofIdx)) {...}
    auto voidToSolidConnections(const std::size_t dofIdx) const
    {
        const auto& ids = voidToSolidConnectionIds().at(dofIdx);
        using Iterator = ConnectionIterator<std::vector<std::size_t>>;
        return Dune::IteratorRange<Iterator>(Iterator(ids.cbegin(), connectionInfo()),
                                             Iterator(ids.cend(), connectionInfo()));
    }

    auto solidToVoidConnections(const std::size_t dofIdx) const
    {
        const auto& ids = solidToVoidConnectionIds().at(dofIdx);
        using Iterator = ConnectionIterator<std::vector<std::size_t>>;
        return Dune::IteratorRange<Iterator>(Iterator(ids.cbegin(), connectionInfo()),
                                             Iterator(ids.cend(), connectionInfo()));
    }

    const auto& voidToSolidConnectionIds() const
    { return voidToSolidConnectionIds_; }

    const auto& solidToVoidConnectionIds() const
    { return solidToVoidConnectionIds_; }

    const auto& connectionInfo() const
    { return connectionInfo_; }

    const auto& voidHostToSubVertexIdxMap() const
    { return voidHostToSubVertexIdxMap_; }

    const auto& solidHostToSubVertexIdxMap() const
    { return solidHostToSubVertexIdxMap_; }

    const auto& voidHostToSubElementIdxMap() const
    { return voidHostToSubElementIdxMap_; }

    const auto& solidHostToSubElementIdxMap() const
    { return solidHostToSubElementIdxMap_; }

    const auto& hostGridElementIndexToGlobalId() const
    { return hostGridElementIndexToGlobalId_; }

private:

    template<class HostGridView, class HostGridData>
    std::vector<HostGridConnectionInfo> getConnectionInfo_(const HostGridView& hostGridView,
                                                           const HostGridData& hostGridData)
    {
        std::vector<HostGridConnectionInfo> connectionInfo;
        std::size_t connectionGlobalId = 0;

        for (const auto& element : elements(hostGridView))
        {
            if (hostGridData.getParameter(element, "ThroatDomainType") == 2) // interconnection throat
            {
                const auto& vertex0 = element.template subEntity<1>(0);
                const auto& vertex1 = element.template subEntity<1>(1);

                HostGridConnectionInfo info;
                info.hostGridElementIndex = hostGridView.indexSet().index(element);
                info.connectionGlobalId = connectionGlobalId++;
                info.voidVertexHostIdx = hostGridView.indexSet().index(vertex0);
                info.solidVertexHostIdx = hostGridView.indexSet().index(vertex1);

                if (hostGridData.getParameter(vertex0, "PoreDomainType") == 1)
                {
                    assert(hostGridData.getParameter(vertex1, "PoreDomainType") == 0);
                    using std::swap;
                    swap(info.voidVertexHostIdx, info.solidVertexHostIdx);
                }

                info.connectionArea = hostGridData.getParameter(element, "ThroatCrossSectionalArea");
                info.connectionLength = element.geometry().volume();

                for (const auto& intersection : intersections(hostGridView, element))
                {
                    if (!intersection.neighbor())
                        continue;

                    const auto& outsideElement = intersection.outside();

                    // skip other interconnection throats
                    if (hostGridData.getParameter(outsideElement, "ThroatDomainType") == 2)
                        continue;

                    const auto outsideElementIdx = hostGridView.indexSet().index(outsideElement);

                    if (hostGridData.getParameter(outsideElement, "ThroatDomainType") == 0)
                        info.voidElementHostIdx.push_back(outsideElementIdx);
                    else
                        info.solidElementHostIdx.push_back(outsideElementIdx);

                    std::array outsideDomainType = {-1, -1};
                    for (int localVIdx = 0; localVIdx < 2; ++localVIdx)
                    {
                        const auto& outsideVertex = outsideElement.template subEntity<1>(localVIdx);
                        const auto outsideVertexIdx = hostGridView.indexSet().index(outsideVertex);
                        outsideDomainType[localVIdx] = hostGridData.getParameter(outsideVertex, "PoreDomainType");

                        if (localVIdx == 1 && (outsideDomainType[1] != outsideDomainType[0]))
                            DUNE_THROW(Dune::IOError, "Pore with hostIdx " << hostGridView.indexSet().index(outsideElement.template subEntity<1>(0))
                                                       << " has domain type " << outsideDomainType[0]
                                                       << ", but pore with hostIdx " << outsideVertexIdx
                                                       << " has domain type " << outsideDomainType[1] << ". Check your grid file");

                        if (outsideDomainType[localVIdx] == 0)
                            info.coupledVoidVertexHostIdx.push_back(outsideVertexIdx);
                        else
                            info.coupledSolidVertexHostIdx.push_back(outsideVertexIdx);
                    }
                }

                connectionInfo.emplace_back(std::move(info));
            }
        }

        return connectionInfo;
    }

    template<class HostGridView, class SubGridGeometry, class Map>
    void fillVertexMap_(const HostGridView& hostGridView, const SubGridGeometry& subGridGeometry, Map& map)
    {
        for (const auto& vertex : vertices(subGridGeometry.gridView()))
        {
            const auto vIdxSub = subGridGeometry.vertexMapper().index(vertex);
            const auto vIdxHost = hostGridView.indexSet().index(vertex.impl().hostEntity());
            map[vIdxHost] = vIdxSub;
        }
    }

    template<class HostGridView, class SubGridGeometry, class Map>
    void fillElementMap_(const HostGridView& hostGridView, const SubGridGeometry& subGridGeometry, Map& map)
    {
        for (const auto& element : elements(subGridGeometry.gridView()))
        {
            const auto eIdxSub = subGridGeometry.elementMapper().index(element);
            const auto eIdxHost = hostGridView.indexSet().index(element.impl().hostEntity());
            map[eIdxHost] = eIdxSub;
        }
    }

    //! Removes duplicate entries from the coupling stencils
    void removeDuplicates_(std::vector<std::size_t>& stencil)
    {
        std::sort(stencil.begin(), stencil.end());
        stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
    }

    std::unordered_map<std::size_t, std::size_t> voidHostToSubVertexIdxMap_;
    std::unordered_map<std::size_t, std::size_t> solidHostToSubVertexIdxMap_;
    std::unordered_map<std::size_t, std::size_t> voidHostToSubElementIdxMap_;
    std::unordered_map<std::size_t, std::size_t> solidHostToSubElementIdxMap_;

    std::vector<bool> isCoupledVoidDof_;
    std::vector<bool> isCoupledSolidDof_;

    std::unordered_map<std::size_t, Stencil> voidToSolidStencils_;
    std::unordered_map<std::size_t, Stencil> solidToVoidStencils_;

    std::unordered_map<std::size_t, std::vector<std::size_t>> voidToSolidConnectionIds_;
    std::unordered_map<std::size_t, std::vector<std::size_t>> solidToVoidConnectionIds_;

    std::vector<SubGridConnectionInfo> connectionInfo_;
    std::unordered_map<std::size_t, std::size_t> hostGridElementIndexToGlobalId_;
};

} // end namespace Dumux::PoreNetwork

#endif
