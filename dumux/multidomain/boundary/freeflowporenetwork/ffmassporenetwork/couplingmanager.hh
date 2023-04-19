// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPoreNetworkCoupling
 * \copydoc Dumux::FreeFlowMassPoreNetworkCouplingManager
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOWPORENETWORK_FFMASSPORENETWORK_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOWPORENETWORK_FFMASSPORENETWORK_COUPLINGMANAGER_HH

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
 * \brief Coupling manager for free-flow mass and pore-network models.
 */
template<class MDTraits>
class FreeFlowMassPoreNetworkCouplingManager
: public CouplingManager<MDTraits>
{
    using Scalar = typename MDTraits::Scalar;
    using ParentType = CouplingManager<MDTraits>;

public:
    static constexpr auto freeFlowMassIndex = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto poreNetworkIndex = typename MDTraits::template SubDomain<1>::Index();

    using SolutionVectorStorage = typename ParentType::SolutionVectorStorage;
private:

    // obtain the type tags of the sub problems
    using FreeFlowMassTypeTag = typename MDTraits::template SubDomain<freeFlowMassIndex>::TypeTag;
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
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename FVElementGeometry<id>::SubControlVolume;

    using VelocityVector = typename Element<freeFlowMassIndex>::Geometry::GlobalCoordinate;

    struct FreeFlowMassCouplingContext
    {
        SubControlVolume<poreNetworkIndex> scv;
        VolumeVariables<poreNetworkIndex> volVars;
        mutable VelocityVector velocity; // velocity needs to be set externally, not available in this class
    };

    struct PoreNetworkCouplingContext
    {
        SubControlVolume<freeFlowMassIndex> scv;
        SubControlVolumeFace<freeFlowMassIndex> scvf;
        VolumeVariables<freeFlowMassIndex> volVars;
        mutable VelocityVector velocity; // velocity needs to be set externally, not available in this class
    };

    using CouplingMapper = StaggeredFreeFlowPoreNetworkCouplingMapper;

public:

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    //! Initialize the coupling manager
    void init(std::shared_ptr<Problem<freeFlowMassIndex>> freeFlowMassProblem,
              std::shared_ptr<Problem<poreNetworkIndex>> pnmProblem,
              std::shared_ptr<CouplingMapper> couplingMapper,
              SolutionVectorStorage& curSol)
    {
        couplingMapper_ = couplingMapper;
        this->setSubProblems(std::make_tuple(freeFlowMassProblem, pnmProblem));
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

        // we need to update all solution-dependent components of the coupling context
        // the dof of domain J has been deflected
        // if domainJ == freeFlowMassIndex: update volvars in the PoreNetworkCouplingContext
        // if domainJ == poreNetworkIndex: update volvars in the FreeFlowMassCouplingContext
        // as the update is symmetric we only need to write this once
        auto& context = std::get<1-j>(couplingContext_);
        for (auto& c : context)
        {
            if (c.scv.dofIndex() == dofIdxGlobalJ)
            {
                const auto& problem = this->problem(domainJ);
                const auto& element = problem.gridGeometry().element(c.scv.elementIndex());
                const auto elemSol = elementSolution(element, this->curSol(domainJ), problem.gridGeometry());
                c.volVars.update(elemSol, problem, element, c.scv);
            }
        }
    }

    // \}

    /*!
     * \brief Access the coupling context needed for the Stokes domain
     */
    const auto& couplingContext(Dune::index_constant<freeFlowMassIndex> domainI,
                                const FVElementGeometry<freeFlowMassIndex>& fvGeometry,
                                const SubControlVolumeFace<freeFlowMassIndex> scvf) const
    {
        auto& contexts = std::get<freeFlowMassIndex>(couplingContext_);
        const auto eIdx = scvf.insideScvIdx();

        if (contexts.empty() || couplingContextBoundForElement_[freeFlowMassIndex] != eIdx)
            bindCouplingContext_(freeFlowMassIndex, fvGeometry);


        return contexts[0];
    }

    /*!
     * \brief Access the coupling context needed for the PNM domain
     */
    const auto& couplingContext(Dune::index_constant<poreNetworkIndex> domainI,
                                const FVElementGeometry<poreNetworkIndex>& fvGeometry,
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
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const Element<i>& element,
                                           Dune::index_constant<j> domainJ) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        if constexpr (domainI == freeFlowMassIndex)
            return couplingMapper_->freeFlowMassToPoreNetworkCouplingStencil(eIdx);
        else
            return couplingMapper_->poreNetworkToFreeFlowMassCouplingStencil(eIdx);
    }

    // \}

    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    bool isCoupled(Dune::index_constant<freeFlowMassIndex> domainI,
                   const SubControlVolumeFace<freeFlowMassIndex>& scvf) const
    { return couplingMapper_->isCoupledFreeFlowMassScvf(scvf.index()); }

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param scv the sub control volume
     */
    bool isCoupled(Dune::index_constant<poreNetworkIndex> domainI,
                   const SubControlVolume<poreNetworkIndex>& scv) const
    { return couplingMapper_->isCoupledPoreNetworkDof(scv.dofIndex()); }

private:
    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    template<std::size_t i>
    bool isCoupledElement_(Dune::index_constant<i> domainI, std::size_t eIdx) const
    {
        if constexpr (domainI == freeFlowMassIndex)
            return couplingMapper_->isCoupledFreeFlowElement(eIdx);
        else
            return couplingMapper_->isCoupledPoreNetworkElement(eIdx);
    }

    //! Return the volume variables of domain i for a given element and scv
    template<std::size_t i>
    VolumeVariables<i> volVars_(Dune::index_constant<i> domainI,
                                const Element<i>& element,
                                const SubControlVolume<i>& scv) const
    {
        VolumeVariables<i> volVars;
        const auto elemSol = elementSolution(
            element, this->curSol(domainI), this->problem(domainI).gridGeometry()
        );
        volVars.update(elemSol, this->problem(domainI), element, scv);
        return volVars;
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
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Darcy element (i.e. Stokes information)
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

        const auto& stencil = couplingStencil(poreNetworkIndex, fvGeometry.element(), freeFlowMassIndex);
        const auto& freeFlowElements = couplingMapper_->pnmElementToFreeFlowElementsMap().at(eIdx);
        auto ffFVGeometry = localView(this->problem(freeFlowMassIndex).gridGeometry());

        for (const auto ffElementIdx : freeFlowElements)
        {
            const auto& ffElement = this->problem(freeFlowMassIndex).gridGeometry().element(ffElementIdx);
            ffFVGeometry.bindElement(ffElement);

            for (const auto& scvf : scvfs(ffFVGeometry))
            {
                if (couplingMapper_->isCoupledFreeFlowMassScvf(scvf.index()))
                {
                    const auto& scv = ffFVGeometry.scv(scvf.insideScvIdx());
                    const auto dofIdx = scv.dofIndex();
                    if (std::any_of(stencil.begin(), stencil.end(), [&](const auto x){ return dofIdx == x; } ))
                    {
                        context.push_back({scv,
                                           scvf,
                                           volVars_(freeFlowMassIndex, ffElement, scv),
                                           VelocityVector{}}
                        );
                    }
                }
            }
        }
    }

        /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Darcy element (i.e. Stokes information)
     */
    void bindCouplingContext_(Dune::index_constant<freeFlowMassIndex> domainI, const FVElementGeometry<freeFlowMassIndex>& fvGeometry) const
    {
        auto& context = std::get<domainI>(couplingContext_);
        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(fvGeometry.element());

        // do nothing if the element is already bound or not coupled to the other domain
        if ((!context.empty() && couplingContextBoundForElement_[domainI] == eIdx) || !isCoupledElement_(domainI, eIdx))
            return;

        context.clear();
        couplingContextBoundForElement_[domainI] = eIdx;


        auto poreNetworkFVGeometry = localView(this->problem(poreNetworkIndex).gridGeometry());
        const auto poreNetworkElemIdx = couplingMapper_->freeFlowElementToPNMElementMap().at(eIdx);
        const auto& poreNetworkElement = this->problem(poreNetworkIndex).gridGeometry().element(poreNetworkElemIdx);
        poreNetworkFVGeometry.bindElement(poreNetworkElement);

        auto poreNetworkScv = [&]
        {
            SubControlVolume<poreNetworkIndex> result;
            std::size_t counter = 0;
            for (auto&& scv : scvs(poreNetworkFVGeometry))
            {
                if (couplingMapper_->isCoupledPoreNetworkDof(scv.dofIndex()))
                {
                    result = scv;
                    ++counter;
                }
            }

            if (counter > 1)
                DUNE_THROW(Dune::InvalidStateException, "Only one pore per throat may be coupled");
            else
                return result;
        }();

        auto volVars = volVars_(poreNetworkIndex, poreNetworkElement, poreNetworkScv);

        context.push_back({std::move(poreNetworkScv),
                           std::move(volVars),
                           VelocityVector{}}
        );
    }

    mutable std::tuple<std::vector<FreeFlowMassCouplingContext>, std::vector<PoreNetworkCouplingContext>> couplingContext_;
    mutable std::array<std::size_t, 2> couplingContextBoundForElement_;

    std::shared_ptr<CouplingMapper> couplingMapper_;
};

} // end namespace Dumux

#endif
