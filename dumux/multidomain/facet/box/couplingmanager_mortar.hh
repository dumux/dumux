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
 * \ingroup FacetCoupling
 * \copydoc Dumux::BoxMortarFacetCouplingManager
 */
#ifndef DUMUX_BOX_MORTAR_FACETCOUPLING_MANAGER_HH
#define DUMUX_BOX_MORTAR_FACETCOUPLING_MANAGER_HH

#include <type_traits>
#include <algorithm>
#include <cassert>
#include <vector>

#include <dune/common/indices.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/promotiontraits.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/evalsolution.hh>

#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/glue.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief TODO
 * \todo TODO: Doc me
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam bulkId The domain id of the bulk problem
 * \tparam facetId The domain id of the lower-dimensional facet problem
 * \tparam mortarId The domain id of the lower-dimensional mortar problem
 */
template< class MDTraits,
          std::size_t bulkId = 0,
          std::size_t facetId = 1,
          std::size_t mortarId = 2>
class BoxMortarFacetCouplingManager
: public virtual CouplingManager< MDTraits >
{
    using ParentType = CouplingManager< MDTraits >;

    using BulkIdType = typename MDTraits::template SubDomain<bulkId>::Index;
    using FacetIdType = typename MDTraits::template SubDomain<facetId>::Index;
    using MortarIdType = typename MDTraits::template SubDomain<mortarId>::Index;

    static constexpr auto bulkDomainId = BulkIdType();
    static constexpr auto facetDomainId = FacetIdType();
    static constexpr auto mortarDomainId = MortarIdType();

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    // further types specific to the sub-problems
    template<std::size_t id> using Scalar = GetPropType<SubDomainTypeTag<id>, Properties::Scalar>;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using NumEqVector = GetPropType<SubDomainTypeTag<id>, Properties::NumEqVector>;
    template<std::size_t id> using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;

    template<std::size_t id> using FVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id> using GridGeometry = FVGridGeometry<id>; // TODO: REMOVE WHEN FIXED
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolume = typename FVGridGeometry<id>::SubControlVolume;
    template<std::size_t id> using SubControlVolumeFace = typename FVGridGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using ElementGeometry = typename Element<id>::Geometry;
    template<std::size_t id> using GlobalPosition = typename ElementGeometry<id>::GlobalCoordinate;
    template<std::size_t id> using GridIndexType = typename GridView<id>::IndexSet::IndexType;
    template<std::size_t id> using Stencil = std::vector<GridIndexType<id>>;

    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using GridVolumeVariables = typename GridVariables<id>::GridVolumeVariables;
    template<std::size_t id> using ElementVolumeVariables = typename GridVolumeVariables<id>::LocalView;

    // promoted coordinate type between two sub-domain grids
    using ctype1 = typename Dune::PromotionTraits< typename GridView<bulkId>::ctype,
                                                   typename GridView<mortarId>::ctype >::PromotedType;
    using ctype = typename Dune::PromotionTraits<typename GridView<facetId>::ctype, ctype1>::PromotedType;

    // extract grid dimensions
    static constexpr int bulkDim = GridView<bulkId>::dimension;
    static constexpr int facetDim = GridView<facetDomainId>::dimension;
    static_assert(facetDim == GridView<mortarDomainId>::dimension, "Grid dimension mismatch");
    static_assert(facetDim == bulkDim-1, "Grid dimension mismatch");

    static constexpr int dimWorld = GridView<bulkDomainId>::dimensionworld;
    static_assert(dimWorld == GridView<facetDomainId>::dimensionworld, "World dimension mismatch");
    static_assert(dimWorld == GridView<mortarDomainId>::dimensionworld, "World dimension mismatch");

    // types to describe the geometry of an interface segment
    using SegmentGeometry = Dune::AffineGeometry<ctype, facetDim, dimWorld>;

    struct BulkMortarInterfaceSegment
    {
        SegmentGeometry geometry;
        GridIndexType<mortarId> mortarElementIndex;
        GridIndexType<bulkId> bulkElementIndex;
        GridIndexType<bulkId> bulkScvfIndex;
        std::size_t sideId;

        BulkMortarInterfaceSegment(const SegmentGeometry& g)
        : geometry(g)
        {}
    };

    struct FacetMortarInterfaceSegment
    {
        SegmentGeometry geometry;
        GridIndexType<mortarId> mortarElementIndex;
        GridIndexType<facetId> facetElementIndex;
        GridIndexType<facetId> facetScvIndex;

        FacetMortarInterfaceSegment(const SegmentGeometry& g)
        : geometry(g)
        {}
    };

    /*!
     * \brief
     * \todo TODO: doc me
     * \note We need unique ptrs here because the local views have no default constructor.
     */
    struct MortarCouplingContext
    {
    public:
        GridIndexType<mortarDomainId> mortarElementIndex;

        // required data in the bulk neighbors
        std::vector< GridIndexType<bulkId> > bulkElementIndices;
        std::vector< std::unique_ptr<FVElementGeometry<bulkId>> > bulkFvGeometries;
        std::vector< std::unique_ptr<ElementVolumeVariables<bulkId>> > bulkElemVolVars;

        // required data in the facet neighbors
        std::vector< GridIndexType<facetId> > facetElementIndices;
        std::vector< std::unique_ptr<FVElementGeometry<facetId>> > facetFvGeometries;
        std::vector< std::unique_ptr<ElementVolumeVariables<facetId>> > facetElemVolVars;

        void clear()
        {
            bulkElementIndices.clear();
            bulkFvGeometries.clear();
            bulkElemVolVars.clear();

            facetElementIndices.clear();
            facetFvGeometries.clear();
            facetElemVolVars.clear();
        }
    };

public:

    //! the type of the solution vector
    using SolutionVector = typename MDTraits::SolutionVector;

    /*!
     * \brief Initialize the coupling manager.
     * \param bulkProblem The problem to be solved on the bulk domain
     * \param facetProblem The problem to be solved on the lower-dimensional domain
     * \param mortarProblem The problem to be solved on the lower-dimensional mortar domain
     * \param curSol The current solution
     */
    void init(std::shared_ptr< Problem<bulkDomainId> > bulkProblem,
              std::shared_ptr< Problem<facetDomainId> > facetProblem,
              std::shared_ptr< Problem<mortarDomainId> > mortarProblem,
              const SolutionVector& curSol)
    {
        // copy the solution vector
        ParentType::updateSolution(curSol);

        this->setSubProblem(bulkProblem, bulkDomainId);
        this->setSubProblem(facetProblem, facetDomainId);
        this->setSubProblem(mortarProblem, mortarDomainId);

        makeBulkCouplingStencils_();
        makeFacetCouplingStencils_();
    }

    /*!
     * \brief The coupling stencil of a given bulk domain element.
     * \note The bulk domain is only coupled to the mortar domain.
     */
    template<std::size_t j>
    const Stencil<j>& couplingStencil(BulkIdType domainId,
                                      const Element<bulkId>& element,
                                      Dune::index_constant<j> domainJ) const
    {
        if constexpr(j == bulkId) { return emptyBulkStencil_; }
        else if constexpr(j == facetId) { return emptyFacetStencil_; }
        else
        {
            static_assert(j == mortarId, "Invalid domain id");

            const auto& bulkGG = this->problem(domainId).fvGridGeometry();
            const auto eIdx = bulkGG.elementMapper().index(element);
            auto it = bulkMortarStencils_.find(eIdx);
            if (it == bulkMortarStencils_.end())
                return emptyMortarStencil_;
            return it->second;
        }
    }

    /*!
     * \brief The coupling stencil of the lower-dimensional facet domain.
     * \note The facet domain is only coupled to the mortar domain.
     */
    template<std::size_t j>
    const Stencil<j>& couplingStencil(FacetIdType domainId,
                                      const Element<facetId>& element,
                                      Dune::index_constant<j> domainJ) const
    {
        if constexpr(j == bulkId) { return emptyBulkStencil_; }
        else if constexpr(j == facetId) { return emptyFacetStencil_; }
        else
        {
            static_assert(j == mortarId, "Invalid domain id");

            const auto& facetGG = this->problem(domainId).fvGridGeometry();
            const auto eIdx = facetGG.elementMapper().index(element);
            return facetMortarStencils_[eIdx];
        }
    }

    /*!
     * \brief The coupling stencil of the mortar domain.
     */
    template<std::size_t j>
    const Stencil<j>& couplingStencil(MortarIdType domainId,
                                      const Element<mortarId>& element,
                                      Dune::index_constant<j> domainJ) const
    {
        const auto& mortarGG = this->problem(domainId).gridGeometry();

        if constexpr(j == bulkId)
        { return mortarBulkStencils_[mortarGG.elementMapper().index(element)]; }
        else if constexpr(j == facetId)
        { return mortarFacetStencils_[mortarGG.elementMapper().index(element)]; }
        else
        {
            static_assert(j == mortarId, "Invalid domain id");
            return emptyMortarStencil_;
        }
    }

    /*!
     * \brief returns true if a bulk scvf lies on a bulk-facet interface.
     * \note for box this is the case if an scvf lies on an interior boundary
     */
    bool isCoupled(const Element<bulkId>& element,
                   const SubControlVolumeFace<bulkId>& scvf) const
    { return scvf.interiorBoundary(); }

    /*!
     * \brief returns true if a bulk scvf lies on a bulk-facet interface.
     */
    bool isOnInteriorBoundary(const Element<bulkId>& element,
                              const SubControlVolumeFace<bulkId>& scvf) const
    { return isCoupled(element, scvf); }

    /*!
     * \brief Evaluates the coupling element residual of a bulk
     *        domain element with respect to a dof of another domain.
     */
    template< class BulkLocalAssembler, std::size_t j >
    typename LocalResidual<bulkId>::ElementResidualVector
    evalCouplingResidual(BulkIdType domainId,
                         const BulkLocalAssembler& bulkLocalAssembler,
                         Dune::index_constant<j> domainIdJ,
                         GridIndexType<j> dofIdxGlobalJ)
    {
        const auto& fvGeometry = bulkLocalAssembler.fvGeometry();
        typename LocalResidual<bulkId>::ElementResidualVector res(fvGeometry.numScv());
        res = 0.0;

        if constexpr(domainIdJ == mortarDomainId)
        {
            const auto& element = bulkLocalAssembler.element();
            const auto eIdx = fvGeometry.fvGridGeometry().elementMapper().index(element);

            auto it = bulkInterfacesMap_.find(eIdx);
            if (it != bulkInterfacesMap_.end())
            {
                for (auto segmentIdx : it->second)
                {
                    const auto& interface = bulkMortarInterfaceSegments_[segmentIdx];
                    const auto mortarElemIdx = interface.mortarElementIndex;
                    const auto& mortarGG = this->problem(mortarDomainId).gridGeometry();
                    const auto mortarElement = mortarGG.element(mortarElemIdx);
                    const auto flux = integrateMortarFlux_(mortarElement,
                                                           mortarElement.geometry(),
                                                           interface.geometry,
                                                           interface.sideId);

                    const auto& scvf = fvGeometry.scvf(interface.bulkScvfIndex);
                    const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
                    res[scv.localDofIndex()] += flux;
                }
            }
        }

        return res;
    }

    /*!
     * \brief Evaluates the coupling element residual of a facet
     *        domain element with respect to a dof of another domain.
     */
    template< class FacetLocalAssembler, std::size_t j>
    typename LocalResidual<facetDomainId>::ElementResidualVector
    evalCouplingResidual(FacetIdType domainId,
                         const FacetLocalAssembler& facetLocalAssembler,
                         Dune::index_constant<j> domainIdJ,
                         GridIndexType<j> dofIdxGlobalJ)
    {
        const auto& fvGeometry = facetLocalAssembler.fvGeometry();
        typename LocalResidual<facetDomainId>::ElementResidualVector res(fvGeometry.numScv());
        res = 0.0;

        if constexpr(domainIdJ == mortarDomainId)
        {
            for (const auto& scv : scvs(facetLocalAssembler.fvGeometry()))
                res[scv.localDofIndex()] -= evalSourcesFromBulk(facetLocalAssembler.element(),
                                                                facetLocalAssembler.fvGeometry(),
                                                                facetLocalAssembler.curElemVolVars(),
                                                                scv);
        }

        return res;
    }

    /*!
     * \brief Evaluates the coupling element residual of a mortar
     *        domain element with respect to a dof of another domain.
     */
    template< class MortarLocalAssembler, std::size_t j >
    typename LocalResidual<mortarDomainId>::ElementResidualVector
    evalCouplingResidual(MortarIdType domainId,
                         const MortarLocalAssembler& mortarLocalAssembler,
                         Dune::index_constant<j> domainIdJ,
                         GridIndexType<j> dofIdxGlobalJ)
    {
        const auto& localResidual = mortarLocalAssembler.localResidual();
        if (GridGeometry<mortarId>::isStandardGalerkin())
            return localResidual.eval(mortarLocalAssembler.element(),
                                      mortarLocalAssembler.feGeometry(),
                                      mortarLocalAssembler.curElemSol());
       return localResidual.eval(mortarLocalAssembler.element(),
                                 mortarLocalAssembler.feGeometry(),
                                 mortarLocalAssembler.curElemSol(),
                                 mortarLocalAssembler.trialSpaceBasisLocalView());
    }

    /*!
     * \brief TODO
     * \todo TODO: Doc me.
     */
    double getMortarFlux(const Element<bulkDomainId>& element,
                         const SubControlVolumeFace<bulkDomainId>& scvf) const
    {
        const auto& bulkGG = this->problem(bulkDomainId).fvGridGeometry();
        const auto eIdx = bulkGG.elementMapper().index(element);

        auto it = bulkInterfacesMap_.find(eIdx);
        if (it != bulkInterfacesMap_.end())
        {
            double flux = 0.0;
            double magnitude = 0.0;
            for (auto interfaceIdx : it->second)
            {
                const auto& interface = bulkMortarInterfaceSegments_[interfaceIdx];
                if (interface.bulkScvfIndex != scvf.index())
                    continue;

                const auto mortarElemIdx = interface.mortarElementIndex;
                const auto& mortarGG = this->problem(mortarDomainId).gridGeometry();
                const auto mortarElement = mortarGG.element(mortarElemIdx);
                flux += integrateMortarFlux_(mortarElement,
                                             mortarElement.geometry(),
                                             interface.geometry,
                                             interface.sideId);
            }

            if (std::abs(magnitude - scvf.area()) > 1e-6)
                DUNE_THROW(Dune::InvalidStateException, "INTEGRATION NOT COMPLET");

            return flux;
        }

        DUNE_THROW(Dune::InvalidStateException, "Could not compute mortar flux");
    }

    /*!
     * \brief Computes the sources (i.e. bulk-facet transfer fluxes)
     *        in a lower-dimensional element of the facet domain.
     */
    NumEqVector<facetDomainId> evalSourcesFromBulk(const Element<facetDomainId>& element,
                                                   const FVElementGeometry<facetDomainId>& fvGeometry,
                                                   const ElementVolumeVariables<facetDomainId>& elemVolVars,
                                                   const SubControlVolume<facetDomainId>& scv)
    {
        NumEqVector<facetDomainId> sources(0.0);

        auto it = facetInterfacesMap_.find(scv.elementIndex());
        if (it == facetInterfacesMap_.end())
            return sources;

        for (auto interfaceIdx : it->second)
        {
            const auto& interface = facetMortarInterfaceSegments_[interfaceIdx];
            if (scv.indexInElement() == interface.facetScvIndex)
            {
                // integrate the mortar variable for each side
                // TODO: Make integration more efficient
                const auto& mortarElemIdx = interface.mortarElementIndex;
                const auto& mortarGG = this->problem(mortarDomainId).gridGeometry();
                const auto mortarElement = mortarGG.element(mortarElemIdx);
                const auto mortarElemGeom = mortarElement.geometry();
                sources[scv.localDofIndex()] -= integrateMortarFlux_(mortarElement, mortarElemGeom,
                                                                     interface.geometry, 0);
                sources[scv.localDofIndex()] -= integrateMortarFlux_(mortarElement, mortarElemGeom,
                                                                     interface.geometry, 1);

            }
        }

        return sources;
    }

    /*!
     * \brief Evaluates the constraint on the mortar at the given position in an element.
     */
    NumEqVector<mortarDomainId> evalMortarCondition(const Element<mortarDomainId>& element,
                                                    const GlobalPosition<mortarDomainId>& globalPos) const
    {
        int hasOnlyNeumann = -1;
        NumEqVector<mortarDomainId> result(0.0);

        const auto& mortarGG = this->problem(mortarDomainId).gridGeometry();
        const auto mortarElemIdx = mortarGG.elementMapper().index(element);

        auto bulkIt = mortarBulkInterfacesMap_.find(mortarElemIdx);
        auto facetIt = mortarFacetInterfacesMap_.find(mortarElemIdx);

        // evaluate fluxes on the bulk sides (TODO: Gravity)
        // on the segment in which the given position lies
        bool addedBulkFlux0 = false;
        bool addedBulkFlux1 = false;
        Scalar<bulkDomainId> interfaceP0 = 0.0;
        Scalar<bulkDomainId> interfaceP1 = 0.0;
        GlobalPosition<bulkDomainId> normal(0.0);

        for (auto interfaceIdx : bulkIt->second)
        {
            const auto& interface = bulkMortarInterfaceSegments_[interfaceIdx];
            if (intersectsPointGeometry(globalPos, interface.geometry))
            {
                if (interface.sideId == 0 && addedBulkFlux0)
                    continue;
                if (interface.sideId == 1 && addedBulkFlux1)
                    continue;

                auto it = std::find(context_.bulkElementIndices.begin(),
                                    context_.bulkElementIndices.end(),
                                    interface.bulkElementIndex);
                assert(it != context_.bulkElementIndices.end());
                const auto i = std::distance(context_.bulkElementIndices.begin(), it);

                const auto& bulkGG = this->problem(bulkDomainId).fvGridGeometry();
                const auto bulkElement = bulkGG.element(interface.bulkElementIndex);
                const auto& fvGeometry = *context_.bulkFvGeometries[i];
                const auto& scvf = fvGeometry.scvf(interface.bulkScvfIndex);
                const auto& bcTypes = this->problem(bulkDomainId).interiorBoundaryTypes(bulkElement, scvf);

                if (hasOnlyNeumann != -1 && bcTypes.hasOnlyNeumann() != bool(hasOnlyNeumann))
                    DUNE_THROW(Dune::InvalidStateException, "Element coupling with heterogeneous conditions!");
                hasOnlyNeumann = int(bcTypes.hasOnlyNeumann());

                const auto elemSol = elementSolution(bulkElement, this->curSol()[bulkDomainId], bulkGG);
                const auto elemGeometry = bulkElement.geometry();
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto& insideVolVars = (*context_.bulkElemVolVars[i])[insideScv];

                // add -n^T*K*gradP and assemble interfaceP
                if (bool(hasOnlyNeumann))
                {
                    const auto gradP = evalGradients(bulkElement, elemGeometry, bulkGG, elemSol, globalPos)[0];
                    result[interface.sideId] -= vtmv(scvf.unitOuterNormal(),
                                                     insideVolVars.permeability(),
                                                     gradP);
                }

                const auto p = evalSolution(bulkElement, elemGeometry, bulkGG, elemSol, globalPos)[0];
                if (interface.sideId == 0) { interfaceP0 = p; addedBulkFlux0 = true; }
                else { normal = scvf.unitOuterNormal(); interfaceP1 = p; addedBulkFlux1 = true; }
            }
        }

        if ( !(addedBulkFlux0 && addedBulkFlux1) )
            DUNE_THROW(Dune::InvalidStateException, "Not all bulk fluxes have been evaluated");

        for (auto interfaceIdx : facetIt->second)
        {
            const auto& interface = facetMortarInterfaceSegments_[interfaceIdx];
            if (intersectsPointGeometry(globalPos, interface.geometry))
            {
                auto it = std::find(context_.facetElementIndices.begin(),
                                    context_.facetElementIndices.end(),
                                    interface.facetElementIndex);
                assert(it != context_.facetElementIndices.end());
                const auto i = std::distance(context_.facetElementIndices.begin(), it);

                const auto& facetGG = this->problem(facetDomainId).fvGridGeometry();
                const auto facetElement = facetGG.element(interface.facetElementIndex);
                const auto elemSol = elementSolution(facetElement, this->curSol()[facetDomainId], facetGG);
                const auto p = evalSolution(facetElement, facetElement.geometry(), facetGG, elemSol, globalPos)[0];

                if (bool(hasOnlyNeumann))
                {
                    const auto& scv = (*context_.facetFvGeometries[i]).scv(interface.facetScvIndex);
                    const auto& volVars = (*context_.facetElemVolVars[i])[scv];
                    const auto l = volVars.extrusionFactor()/2.0;
                    const auto nKn = vtmv(normal, volVars.permeability(), normal);
                    result[0] -= nKn*(p - interfaceP0)/l;
                    result[1] -= nKn*(p - interfaceP1)/l;
                }
                else
                {
                    result[0] = interfaceP0 - p;
                    result[1] = interfaceP1 - p;
                }
            }
        }

        return result;
    }

    using ParentType::bindCouplingContext;
    /*!
     * \brief TODO
     * \todo TODO: Doc me.
     */
    template< class Assembler >
    void bindCouplingContext(MortarIdType domainId,
                             const Element<mortarDomainId>& element,
                             const Assembler& assembler)
    {
        context_.clear();

        const auto& mortarGG = this->problem(mortarDomainId).gridGeometry();
        const auto mortarElemIdx = mortarGG.elementMapper().index(element);

        auto bulkIt = mortarBulkInterfacesMap_.find(mortarElemIdx);
        auto facetIt = mortarFacetInterfacesMap_.find(mortarElemIdx);

        if (bulkIt == mortarBulkInterfacesMap_.end()
            && facetIt == mortarFacetInterfacesMap_.end())
            return;

        // bulk data preparation
        for (auto interfaceIdx : bulkIt->second)
        {
            const auto& interface = bulkMortarInterfaceSegments_[interfaceIdx];
            if (std::count(context_.bulkElementIndices.begin(),
                           context_.bulkElementIndices.end(),
                           interface.bulkElementIndex))
                continue;

            const auto& bulkGG = this->problem(bulkDomainId).fvGridGeometry();
            const auto bulkElement = bulkGG.element(interface.bulkElementIndex);

            auto fvGeometry = localView(bulkGG);
            auto elemVolVars = localView(assembler.gridVariables(bulkDomainId).curGridVolVars());

            fvGeometry.bind(bulkElement);
            elemVolVars.bind(bulkElement, fvGeometry, this->curSol()[bulkDomainId]);

            context_.bulkFvGeometries.emplace_back(std::make_unique<FVElementGeometry<bulkId>>(std::move(fvGeometry)) );
            context_.bulkElemVolVars.emplace_back(std::make_unique<ElementVolumeVariables<bulkId>>(std::move(elemVolVars)) );
        }

        // facet data preparation
        for (auto interfaceIdx : facetIt->second)
        {
            const auto& interface = facetMortarInterfaceSegments_[interfaceIdx];
            if (std::count(context_.facetElementIndices.begin(),
                           context_.facetElementIndices.end(),
                           interface.facetElementIndex))
                continue;

            const auto& facetGG = this->problem(facetDomainId).fvGridGeometry();
            const auto facetElement = facetGG.element(interface.facetElementIndex);

            auto fvGeometry = localView(facetGG);
            auto elemVolVars = localView(assembler.gridVariables(facetDomainId).curGridVolVars());

            fvGeometry.bind(facetElement);
            elemVolVars.bind(facetElement, fvGeometry, this->curSol()[facetDomainId]);

            context_.facetFvGeometries.emplace_back(std::make_unique<FVElementGeometry<facetId>>(std::move(fvGeometry)) );
            context_.facetElemVolVars.emplace_back(std::make_unique<ElementVolumeVariables<facetId>>(std::move(elemVolVars)) );
        }
    }

    /*!
     * \brief TODO DOC THIS OVERLOAD
     * \todo TODO: instead of a complete re-bind, be more selective/efficient!
     */
    template< std::size_t i, class LocalAssembler, std::size_t j,
              std::enable_if_t<i != mortarId, int> = 0>
    void updateCouplingContext(Dune::index_constant<i> domainId,
                               const LocalAssembler& mortarLocalAssembler,
                               Dune::index_constant<j> domainIdJ,
                               GridIndexType<j> dofIdxGlobalJ,
                               const PrimaryVariables<j>& priVarsJ,
                               unsigned int pvIdxJ)
    {}

    /*!
     * \brief After deflecting the solution of a sub-domain during
     *        assembly of the mortar domain, update the context.
     * \todo TODO: instead of a complete re-bind, be more selective/efficient!
     */
    template< class MortarLocalAssembler, std::size_t j >
    void updateCouplingContext(MortarIdType domainId,
                               const MortarLocalAssembler& mortarLocalAssembler,
                               Dune::index_constant<j> domainIdJ,
                               GridIndexType<j> dofIdxGlobalJ,
                               const PrimaryVariables<j>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate deflected solution
        ParentType::updateCouplingContext(domainId, mortarLocalAssembler, domainIdJ,
                                          dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // update the context
        const auto& bulkGG = this->problem(bulkDomainId).fvGridGeometry();
        for (unsigned int i = 0; i < context_.bulkElementIndices.size(); ++i)
        {
            const auto bulkElement = bulkGG.element(context_.bulkElementIndices[i]);
            const auto& bulkFvGeometry = *context_.bulkFvGeometries[i];
            (*context_.bulkElemVolVars[i]).bind(bulkElement, bulkFvGeometry, this->curSol()[bulkDomainId]);
        }

        const auto& facetGG = this->problem(facetDomainId).fvGridGeometry();
        for (unsigned int i = 0; i < context_.facetElementIndices.size(); ++i)
        {
            const auto facetElement = facetGG.element(context_.facetElementIndices[i]);
            const auto& facetFvGeometry = *context_.facetFvGeometries[i];
            (*context_.facetElemVolVars[i]).bind(facetElement, facetFvGeometry, this->curSol()[facetDomainId]);
        }
    }

private:
    //! Intersect bulk and mortar grids and determine coupling stencils
    void makeBulkCouplingStencils_()
    {
        const auto& bulkGG = this->problem(bulkDomainId).fvGridGeometry();
        const auto& mortarGG = this->problem(mortarDomainId).gridGeometry();
        const auto glue = makeGlue(mortarGG, bulkGG);

        // determine estimate for number of interfaces
        std::size_t numInterfaces = 0;
        for (const auto& is : intersections(glue))
            numInterfaces += is.numTargetNeighbors();

        bulkMortarStencils_.clear();
        mortarBulkStencils_.clear();
        bulkMortarInterfaceSegments_.clear();
        mortarBulkStencils_.resize(mortarGG.gridView().size(0));
        bulkMortarInterfaceSegments_.reserve(numInterfaces);

        for (const auto& is : intersections(glue))
        {
            if (is.numDomainNeighbors() != 1)
                DUNE_THROW(Dune::InvalidStateException, "More than one mortar element in intersection");

            const auto mortarElement = is.domainEntity(0);
            const auto mortarElementIdx = mortarGG.elementMapper().index(mortarElement);
            const auto mortarElemGeometry = mortarElement.geometry();

            auto mortarFeGeometry = localView(mortarGG);
            mortarFeGeometry.bind(mortarElement);

            std::vector< GridIndexType<mortarDomainId> > mortarElemDofs;
            for (unsigned int i = 0; i < mortarFeGeometry.feBasisLocalView().size(); ++i)
                mortarElemDofs.push_back(mortarFeGeometry.feBasisLocalView().index(i));

            for (unsigned int nIdx = 0; nIdx < is.numTargetNeighbors(); ++nIdx)
            {
                const auto bulkElement = is.targetEntity(nIdx);
                const auto bulkElementIdx = bulkGG.elementMapper().index(bulkElement);

                auto fvGeometry = localView(bulkGG);
                fvGeometry.bindElement(bulkElement);

                // fill stencils
                for (const auto& scv : scvs(fvGeometry))
                    mortarBulkStencils_[mortarElementIdx].push_back(scv.dofIndex());
                for (auto mortarDofIdx : mortarElemDofs)
                    bulkMortarStencils_[bulkElementIdx].push_back(mortarDofIdx);

                for (const auto& scvf : scvfs(fvGeometry))
                {
                    using ScvfGeometry = typename SubControlVolumeFace<bulkId>::Traits::Geometry;
                    using MortarElemGeometry = typename Element<mortarId>::Geometry;
                    using ISAlgorithm = GeometryIntersection<ScvfGeometry, MortarElemGeometry>;

                    typename ISAlgorithm::Intersection result;
                    if (ISAlgorithm::intersection(scvf.geometry(), mortarElemGeometry, result))
                    {
                        const auto geometries = makeInterfaceSegmentGeometries_(result);
                        for (const auto& geometry : geometries)
                        {
                            // create new interface
                            const auto iFaceIdx = bulkMortarInterfaceSegments_.size();
                            bulkMortarInterfaceSegments_.emplace_back(geometry);

                            auto& interface = bulkMortarInterfaceSegments_.back();
                            interface.bulkElementIndex = bulkElementIdx;
                            interface.bulkScvfIndex = scvf.index();
                            interface.mortarElementIndex = mortarElementIdx;
                            interface.sideId = nIdx;

                            bulkInterfacesMap_[bulkElementIdx].push_back(iFaceIdx);
                            mortarBulkInterfacesMap_[mortarElementIdx].push_back(iFaceIdx);
                        }
                    }
                }
            }
        }

        auto makeUnique = [] (auto& s)
        {
            std::sort(s.begin(), s.end());
            s.erase(std::unique(s.begin(), s.end()), s.end());
        };
        for (auto& dataPair : bulkMortarStencils_) makeUnique(dataPair.second);
        for (auto& s : mortarBulkStencils_) makeUnique(s);
    }

    //! Intersect facet and mortar grids and determine coupling stencils
    void makeFacetCouplingStencils_()
    {
        const auto& facetGG = this->problem(facetDomainId).fvGridGeometry();
        const auto& mortarGG = this->problem(mortarDomainId).gridGeometry();
        const auto glue = makeGlue(mortarGG, facetGG);

        // determine estimate for number of interfaces
        std::size_t numInterfaces = 0;
        for (const auto& is : intersections(glue))
            numInterfaces += is.numTargetNeighbors();

        facetMortarStencils_.clear();
        mortarFacetStencils_.clear();
        facetMortarInterfaceSegments_.clear();
        facetMortarStencils_.resize(facetGG.gridView().size(0));
        facetMortarInterfaceSegments_.reserve(numInterfaces);

        for (const auto& is : intersections(glue))
        {
            if (is.numDomainNeighbors() != 1)
                DUNE_THROW(Dune::InvalidStateException, "More than one mortar element in intersection");
            if (is.numTargetNeighbors() != 1)
                DUNE_THROW(Dune::InvalidStateException, "More than one facet element in intersection");

            const auto mortarElement = is.domainEntity(0);
            const auto facetElement = is.targetEntity(0);
            const auto mortarElementIdx = mortarGG.elementMapper().index(mortarElement);
            const auto facetElementIdx = facetGG.elementMapper().index(facetElement);
            const auto mortarElemGeometry = mortarElement.geometry();

            auto mortarFeGeometry = localView(mortarGG);
            auto facetFvGeometry = localView(facetGG);
            mortarFeGeometry.bind(mortarElement);
            facetFvGeometry.bindElement(facetElement);

            for (unsigned int i = 0; i < mortarFeGeometry.feBasisLocalView().size(); ++i)
                facetMortarStencils_[facetElementIdx].push_back(mortarFeGeometry.feBasisLocalView().index(i));

            for (const auto& scv : scvs(facetFvGeometry))
            {
                mortarFacetStencils_[mortarElementIdx].push_back(scv.dofIndex());

                using ScvGeometry = typename SubControlVolume<facetId>::Traits::Geometry;
                using MortarElemGeometry = typename Element<mortarId>::Geometry;
                using ISAlgorithm = GeometryIntersection<ScvGeometry, MortarElemGeometry>;

                typename ISAlgorithm::Intersection result;
                if (ISAlgorithm::intersection(scv.geometry(), mortarElemGeometry, result))
                {
                    const auto geometries = makeInterfaceSegmentGeometries_(result);
                    for (const auto& geometry : geometries)
                    {
                        const auto iFaceIdx = facetMortarInterfaceSegments_.size();
                        facetMortarInterfaceSegments_.emplace_back(geometry);

                        auto& interface = facetMortarInterfaceSegments_.back();
                        interface.facetElementIndex = facetElementIdx;
                        interface.facetScvIndex = scv.indexInElement();
                        interface.mortarElementIndex = mortarElementIdx;

                        facetInterfacesMap_[facetElementIdx].push_back(iFaceIdx);
                        mortarFacetInterfacesMap_[mortarElementIdx].push_back(iFaceIdx);
                    }
                }
            }
        }

        auto makeUnique = [] (auto& stencil)
        {
            std::sort(stencil.begin(), stencil.end());
            stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
        };
        for (auto& s : mortarFacetStencils_) makeUnique(s);
        for (auto& s : facetMortarStencils_) makeUnique(s);
    }

    //! turn a 1d raw entity intersection into a segment geometry
    template<class RawIntersection, int d = facetDim, std::enable_if_t<d == 1, int> = 0>
    std::array<SegmentGeometry, 1> makeInterfaceSegmentGeometries_(const RawIntersection& is) const
    { return {SegmentGeometry(Dune::GeometryTypes::line, is)}; }

    //! turn a 2d raw entity intersection into a segment geometry
    template<class RawIntersection, int d = facetDim, std::enable_if_t<d == 2, int> = 0>
    std::vector<SegmentGeometry> makeInterfaceSegmentGeometries_(const RawIntersection& is) const
    {
        std::vector<SegmentGeometry> result;
        for (const auto& triangle : triangulate<facetDim, dimWorld>(is))
            result.emplace_back(Dune::GeometryTypes::triangle, triangle);
        return result;
    }

    //! TODO: Doc me
    Scalar<mortarDomainId> integrateMortarFlux_(const Element<mortarDomainId>& element,
                                                const ElementGeometry<mortarDomainId>& geometry,
                                                const SegmentGeometry& segmentGeometry,
                                                unsigned int pvIdx) const
    {
        const auto& mortarGG = this->problem(mortarDomainId).gridGeometry();
        const auto& mortarSol = this->curSol()[mortarDomainId];
        const auto elemSol = elementSolution(element, mortarSol, mortarGG);

        Scalar<mortarDomainId> flux = 0.0;
        static const auto intOrder = getParamFromGroup<Scalar<mortarDomainId>>("Mortar", "IntegrationOrder");

        using QuadRules = Dune::QuadratureRules<Scalar<mortarDomainId>, facetDim>;
        for (const auto& qp : QuadRules::rule(segmentGeometry.type(), intOrder))
        {
            const auto ipLocal = qp.position();
            const auto ipGlobal = geometry.global(ipLocal);
            const auto ipSol = evalSolution(element, geometry, mortarGG, elemSol, ipGlobal);
            flux += ipSol[pvIdx]*qp.weight()*geometry.integrationElement(ipLocal);
        }

        return flux;
    }

    //! interfaces
    std::vector<BulkMortarInterfaceSegment> bulkMortarInterfaceSegments_;
    std::vector<FacetMortarInterfaceSegment> facetMortarInterfaceSegments_;

    //! coupling stencils
    std::unordered_map< GridIndexType<bulkId>, Stencil<mortarId> > bulkMortarStencils_;
    std::vector< Stencil<mortarId> > facetMortarStencils_;

    std::vector< Stencil<bulkId> > mortarBulkStencils_;
    std::vector< Stencil<facetId> > mortarFacetStencils_;

    //! Maps to the elements the embedded interfaces
    using IdxList = std::vector<std::size_t>;
    std::unordered_map< GridIndexType<bulkId>, IdxList > bulkInterfacesMap_;
    std::unordered_map< GridIndexType<facetId>, IdxList > facetInterfacesMap_;
    std::unordered_map< GridIndexType<mortarId>, IdxList > mortarBulkInterfacesMap_;
    std::unordered_map< GridIndexType<mortarId>, IdxList > mortarFacetInterfacesMap_;

    //! empty stencils to return for non-coupled elements
    Stencil<bulkId> emptyBulkStencil_;
    Stencil<facetId> emptyFacetStencil_;
    Stencil<mortarId> emptyMortarStencil_;

    //! The coupling contexts for the mortar domain
    MortarCouplingContext context_;
};

} // end namespace Dumux

#endif // DUMUX_BOX_MORTAR_FACETCOUPLING_MANAGER_HH
