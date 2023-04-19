// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPoreNetworkCoupling
 * \copydoc Dumux::FreeFlowMomentumPoreNetworkCouplingManager
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMMOMENTUMPORENETWORK_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMMOMENTUMPORENETWORK_COUPLINGMANAGER_HH

#include <utility>
#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/boundary/freeflowporenetwork/couplingmapper.hh>

namespace Dumux {

/*!
 * \ingroup FreeFlowPoreNetworkCoupling
 * \brief Coupling manager for free-flow momentum and pore-network models.
 */
template<class MDTraits>
class FreeFlowMomentumPoreNetworkCouplingManager
: public CouplingManager<MDTraits>
{
    using Scalar = typename MDTraits::Scalar;
    using ParentType = CouplingManager<MDTraits>;

public:
    static constexpr auto freeFlowMomentumIndex = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto poreNetworkIndex = typename MDTraits::template SubDomain<1>::Index();

    using SolutionVectorStorage = typename ParentType::SolutionVectorStorage;
private:
    // obtain the type tags of the sub problems
    using FreeFlowMomentumTypeTag = typename MDTraits::template SubDomain<freeFlowMomentumIndex>::TypeTag;
    using PoreNetworkTypeTag = typename MDTraits::template SubDomain<poreNetworkIndex>::TypeTag;

    using CouplingStencils = std::unordered_map<std::size_t, std::vector<std::size_t> >;
    using CouplingStencil = CouplingStencils::mapped_type;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using GridView = typename GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>::GridView;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using GridVolumeVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>;
    template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using GridFluxVariablesCache = GetPropType<SubDomainTypeTag<id>, Properties::GridFluxVariablesCache>;
    template<std::size_t id> using ElementFluxVariablesCache = typename GridFluxVariablesCache<id>::LocalView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename FVElementGeometry<id>::SubControlVolume;

    using VelocityVector = typename Element<freeFlowMomentumIndex>::Geometry::GlobalCoordinate;

    struct FreeFlowMomentumCouplingContext
    {
        FVElementGeometry<poreNetworkIndex> fvGeometry;
        ElementVolumeVariables<poreNetworkIndex> elemVolVars;
        ElementFluxVariablesCache<poreNetworkIndex> elemFluxVarsCache;
        std::size_t poreNetworkDofIdx;
    };

    struct PoreNetworkCouplingContext
    {
        SubControlVolumeFace<freeFlowMomentumIndex> freeFlowMomentumScvf;
        VelocityVector faceVelocity;
        std::size_t freeFlowMomentumDofIdx;
    };

    using CouplingMapper = StaggeredFreeFlowPoreNetworkCouplingMapper;
    using GridVariablesTuple = typename MDTraits::template TupleOfSharedPtr<GridVariables>;

public:

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    //! Initialize the coupling manager
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> freeFlowMomentumProblem,
              std::shared_ptr<Problem<poreNetworkIndex>> porousMediumProblem,
              GridVariablesTuple&& gridVariables,
              std::shared_ptr<CouplingMapper> couplingMapper,
              SolutionVectorStorage& curSol)
    {
        couplingMapper_ = couplingMapper;
        gridVariables_ = gridVariables;
        this->setSubProblems(std::make_tuple(freeFlowMomentumProblem, porousMediumProblem));
        this->attachSolution(curSol);
    }

    // \}


    /*!
     * \brief Methods to be accessed by the assembly
     */
    // \{

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual (called from the local assembler)
     */
    template<std::size_t i, class Assembler>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<i>& element, const Assembler& assembler) const
    {
        bindCouplingContext_(domainI, element);
    }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual (called from the local assembler)
     */
    template<std::size_t i>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<i>& element) const
    {
        bindCouplingContext_(domainI, element);
    }

    /*!
     * \brief Update the coupling context
     */
    template<std::size_t i, std::size_t j, class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<j> domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<j>& priVarsJ,
                               int pvIdxJ)
    {
        this->curSol(domainJ)[dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        const auto eIdx = localAssemblerI.fvGeometry().gridGeometry().elementMapper().index(localAssemblerI.fvGeometry().element());
        if (!isCoupledElement_(domainI, eIdx))
            return;

        // we need to update all solution-depenent components of the coupling context
        // the dof of domain J has been deflected

        // update the faceVelocity in the PoreNetworkCouplingContext
        if constexpr (domainJ == freeFlowMomentumIndex)
        {
            // we only need to update if we are assembling the porous medium domain
            // since the freeflow domain will not use the velocity from the context
            if constexpr (domainI == poreNetworkIndex)
            {
                auto& context = std::get<poreNetworkIndex>(couplingContext_);
                for (auto& c : context)
                {
                    if (c.freeFlowMomentumDofIdx == dofIdxGlobalJ)
                    {
                        assert(c.freeFlowMomentumScvf.isFrontal() && c.freeFlowMomentumScvf.boundary());
                        c.faceVelocity = faceVelocity(c.freeFlowMomentumScvf, c.freeFlowMomentumDofIdx);
                    }
                }
            }
        }

        // update the elemVolVars and elemFluxVarsCache in the FreeFlowMomentumCouplingContext
        else if constexpr (domainJ == poreNetworkIndex)
        {
            assert(couplingContextBoundForElement_[domainI] == localAssemblerI.fvGeometry().gridGeometry().elementMapper().index(localAssemblerI.fvGeometry().element()));
            // there is only one context per coupled free-flow momentum dof
            auto& context = std::get<freeFlowMomentumIndex>(couplingContext_)[0];
            const auto& ggJ = context.fvGeometry.gridGeometry();
            const auto& element = context.fvGeometry.element();
            const auto elemSol = elementSolution(element, this->curSol(domainJ), ggJ);

            for (const auto& scv : scvs(context.fvGeometry))
            {
                if (scv.dofIndex() == dofIdxGlobalJ)
                {
                    if constexpr (ElementVolumeVariables<poreNetworkIndex>::GridVolumeVariables::cachingEnabled)
                        gridVars_(poreNetworkIndex).curGridVolVars().volVars(scv).update(std::move(elemSol), this->problem(domainJ), element, scv);
                    else
                        context.elemVolVars[scv].update(std::move(elemSol), this->problem(domainJ), element, scv);
                }
            }

            const auto& scvf = context.fvGeometry.scvf(0);
            if constexpr (ElementFluxVariablesCache<poreNetworkIndex>::GridFluxVariablesCache::cachingEnabled)
            {
                const auto eIdx = ggJ.elementMapper().index(element);
                gridVars_(poreNetworkIndex).gridFluxVarsCache().cache(eIdx, scvf.index()).update(this->problem(domainJ), element, context.fvGeometry, context.elemVolVars, scvf);
            }
            else
                context.elemFluxVarsCache[scvf].update(this->problem(domainJ), element, context.fvGeometry, context.elemVolVars, scvf);
        }
    }

    // \}

    /*!
     * \brief Access the coupling context needed for the Stokes domain
     */
    const auto& couplingContext(const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                const SubControlVolumeFace<freeFlowMomentumIndex> scvf) const
    {
        auto& contexts = std::get<freeFlowMomentumIndex>(couplingContext_);

        if (contexts.empty() || couplingContextBoundForElement_[freeFlowMomentumIndex] != fvGeometry.elementIndex())
            bindCouplingContext_(freeFlowMomentumIndex, fvGeometry);


        return contexts[0];
    }

    /*!
     * \brief Access the coupling context needed for the PNM domain
     */
    const auto& couplingContext(const FVElementGeometry<poreNetworkIndex>& fvGeometry,
                                const SubControlVolume<poreNetworkIndex> scv) const
    {
        auto& contexts = std::get<poreNetworkIndex>(couplingContext_);

        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(fvGeometry.element());
        if (contexts.empty() || couplingContextBoundForElement_[poreNetworkIndex] != eIdx)
            bindCouplingContext_(poreNetworkIndex, fvGeometry);

        return contexts;
    }

    /*!
     * \brief The coupling stencils
     */
    // \{

    /*!
     * \brief The Stokes cell center coupling stencil w.r.t. Darcy DOFs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<poreNetworkIndex> domainI,
                                           const Element<poreNetworkIndex>& element,
                                           Dune::index_constant<freeFlowMomentumIndex> domainJ) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        return couplingMapper_->poreNetworkToFreeFlowMomentumCouplingStencil(eIdx);
    }

    /*!
     * \brief returns an iterable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the residual of the given sub-control volume of domain i
     *
     * \param domainI the domain index of domain i
     * \param elementI the coupled element of domain í
     * \param scvI the sub-control volume of domain i
     * \param domainJ the domain index of domain j
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                           const Element<freeFlowMomentumIndex>& elementI,
                                           const SubControlVolume<freeFlowMomentumIndex>& scvI,
                                           Dune::index_constant<poreNetworkIndex> domainJ) const
    {
        return couplingMapper_->freeFlowMomentumToPoreNetworkCouplingStencil(scvI.dofIndex());
    }

    using ParentType::evalCouplingResidual;

    /*!
     * \brief evaluate the coupling residual
     * special interface for fcstaggered methods
     */
    template<class LocalAssemblerI, std::size_t j>
    decltype(auto) evalCouplingResidual(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        const SubControlVolume<freeFlowMomentumIndex>& scvI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ) const
    {
        return localAssemblerI.evalLocalResidual();
    }

    // \}

    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    bool isCoupledLateralScvf(Dune::index_constant<freeFlowMomentumIndex> domainI, const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    { return  couplingMapper_->isCoupledFreeFlowMomentumLateralScvf(scvf.index()); }

    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    bool isCoupled(Dune::index_constant<freeFlowMomentumIndex> domainI, const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        return couplingMapper_->isCoupledFreeFlowMomentumScvf(scvf.index()) || couplingMapper_->isCoupledFreeFlowMomentumLateralScvf(scvf.index());
    }

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param scv the sub control volume
     */
    bool isCoupled(Dune::index_constant<poreNetworkIndex> domainI,
                   const SubControlVolume<poreNetworkIndex>& scv) const
    { return couplingMapper_->isCoupledPoreNetworkDof(scv.dofIndex()); }

    /*!
     * \brief Returns the velocity at a given sub control volume face.
     */
    auto faceVelocity(const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                      std::size_t freeFlowMomentumDofIdx) const
    {
        // create a unit normal vector oriented in positive coordinate direction
        auto velocity = scvf.unitOuterNormal();
        using std::abs;
        std::for_each(velocity.begin(), velocity.end(), [](auto& v){ v = abs(v); });

        // create the actual velocity vector
        velocity *= this->curSol(freeFlowMomentumIndex)[freeFlowMomentumDofIdx];

        return velocity;
    }

private:

    /*!
     * \brief Returns whether a given element is coupled to the other domain
     */
    template<std::size_t i>
    bool isCoupledElement_(Dune::index_constant<i> domainI, std::size_t eIdx) const
    {
        if constexpr (i == freeFlowMomentumIndex)
            return couplingMapper_->isCoupledFreeFlowElement(eIdx);
        else
            return couplingMapper_->isCoupledPoreNetworkElement(eIdx);
    }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Darcy element (i.e. Stokes information)
     */
    template<std::size_t i>
    void bindCouplingContext_(Dune::index_constant<i> domainI, const Element<i>& element) const
    {
        const auto fvGeometry = localView(this->problem(domainI).gridGeometry()).bindElement(element);
        bindCouplingContext_(domainI, fvGeometry);
    }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an free-flow momentum element
     */
    void bindCouplingContext_(Dune::index_constant<freeFlowMomentumIndex> domainI, const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry) const
    {
        auto& context = std::get<domainI>(couplingContext_);
        const auto eIdx = fvGeometry.elementIndex();

        // do nothing if the element is already correctly bound
        if ((!context.empty() && couplingContextBoundForElement_[domainI] == eIdx))
            return;

        bool bindElement = false;
        std::size_t actuallyCoupledFreeFlowElementIndex;

        // if the element is directly coupled to a lowDim dof,
        // bind the element itself
        if (couplingMapper_->isCoupledFreeFlowElement(eIdx))
        {
            bindElement = true;
            actuallyCoupledFreeFlowElementIndex = eIdx;
        }
        else
        {
            // if we assemble another element that is not directly coupled to the lowDim dof
            // but shares an intersection (and hence, a dof) with a neighbor element that does, bind that neighbor
            for (const auto& intersection : intersections(fvGeometry.gridGeometry().gridView(), fvGeometry.element()))
            {
                const auto dofIdx = fvGeometry.gridGeometry().intersectionMapper().globalIntersectionIndex(fvGeometry.element(), intersection.indexInInside());
                if (couplingMapper_->isCoupledFreeFlowMomentumDof(dofIdx))
                {
                    bindElement = true;
                    actuallyCoupledFreeFlowElementIndex = fvGeometry.gridGeometry().elementMapper().index(intersection.outside());
                }
            }
        }

        // do nothing if the element is not coupled to the other domain
        if (!bindElement)
            return;

        context.clear();
        couplingContextBoundForElement_[domainI] = eIdx;

        auto poreNetworkFVGeometry = localView(this->problem(poreNetworkIndex).gridGeometry());
        auto poreNetworkElemVolVars = localView(gridVars_(poreNetworkIndex).curGridVolVars());
        auto poreNetworkElemFluxVarsCache = localView(gridVars_(poreNetworkIndex).gridFluxVarsCache());

        const auto poreNetworkElemIdx = couplingMapper_->freeFlowElementToPNMElementMap().at(actuallyCoupledFreeFlowElementIndex);
        const auto& poreNetworkElement = this->problem(poreNetworkIndex).gridGeometry().element(poreNetworkElemIdx);

        poreNetworkFVGeometry.bindElement(poreNetworkElement);
        poreNetworkElemVolVars.bind(poreNetworkElement, poreNetworkFVGeometry, this->curSol(poreNetworkIndex));
        poreNetworkElemFluxVarsCache.bind(poreNetworkElement, poreNetworkFVGeometry, poreNetworkElemVolVars);

        const std::size_t poreNetworkDofIdx = [&]
        {
            std::size_t idx = 0;
            std::size_t counter = 0;
            for (const auto& scv : scvs(poreNetworkFVGeometry))
            {
                if (couplingMapper_->isCoupledPoreNetworkDof(scv.dofIndex()))
                {
                    idx = scv.dofIndex();
                    ++counter;
                }
            }

            if (counter != 1)
                DUNE_THROW(Dune::InvalidStateException, "Exactly one pore per throat needs to be coupled with the FF domain");

            return idx;
        }();

        context.push_back({std::move(poreNetworkFVGeometry),
                           std::move(poreNetworkElemVolVars),
                           std::move(poreNetworkElemFluxVarsCache),
                           poreNetworkDofIdx}
        );
    }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an pore-network model element
     */
    void bindCouplingContext_(Dune::index_constant<poreNetworkIndex> domainI, const FVElementGeometry<poreNetworkIndex>& fvGeometry) const
    {
        auto& context = std::get<domainI>(couplingContext_);
        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(fvGeometry.element());

        // do nothing if the element is already bound or not coupled to the other domain
        if ((!context.empty() && couplingContextBoundForElement_[domainI] == eIdx) || !isCoupledElement_(domainI, eIdx))
            return;

        context.clear();
        couplingContextBoundForElement_[domainI] = eIdx;

        const auto& stencil = couplingStencil(poreNetworkIndex, fvGeometry.element(), freeFlowMomentumIndex);
        const auto& freeFlowElements = couplingMapper_->pnmElementToFreeFlowElementsMap().at(eIdx);
        auto ffFVGeometry = localView(this->problem(freeFlowMomentumIndex).gridGeometry());

        for (const auto ffElementIdx : freeFlowElements)
        {
            const auto& ffElement = this->problem(freeFlowMomentumIndex).gridGeometry().element(ffElementIdx);
            ffFVGeometry.bindElement(ffElement);
            for (const auto& scv : scvs(ffFVGeometry))
            {
                if (couplingMapper_->isCoupledFreeFlowMomentumDof(scv.dofIndex()))
                {
                    if (std::any_of(stencil.begin(), stencil.end(), [&](const auto x){ return scv.dofIndex() == x; } ))
                    {
                        const auto& coupledScvf = ffFVGeometry.frontalScvfOnBoundary(scv);
                        context.push_back({coupledScvf,
                                           faceVelocity(coupledScvf, scv.dofIndex()),
                                           scv.dofIndex()}
                        );
                    }
                }
            }
        }
    }

    /*!
     * \brief Return a reference to the grid variables of a sub problem
     * \param domainIdx The domain index
     */
    template<std::size_t i>
    const GridVariables<i>& gridVars_(Dune::index_constant<i> domainIdx) const
    {
        if (std::get<i>(gridVariables_))
            return *std::get<i>(gridVariables_);
        else
            DUNE_THROW(Dune::InvalidStateException, "The gridVariables pointer was not set. Use setGridVariables() before calling this function");
    }

    /*!
     * \brief Return a reference to the grid variables of a sub problem
     * \param domainIdx The domain index
     */
    template<std::size_t i>
    GridVariables<i>& gridVars_(Dune::index_constant<i> domainIdx)
    {
        if (std::get<i>(gridVariables_))
            return *std::get<i>(gridVariables_);
        else
            DUNE_THROW(Dune::InvalidStateException, "The gridVariables pointer was not set. Use setGridVariables() before calling this function");
    }

    //! A tuple of std::shared_ptrs to the grid variables of the sub problems
    GridVariablesTuple gridVariables_;

    mutable std::tuple<std::vector<FreeFlowMomentumCouplingContext>, std::vector<PoreNetworkCouplingContext>> couplingContext_;
    mutable std::array<std::size_t, 2> couplingContextBoundForElement_;

    std::shared_ptr<CouplingMapper> couplingMapper_;
};

} // end namespace Dumux

#endif
