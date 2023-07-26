// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FFMassPMCouplingManagerStaggeredCCTpfa
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMASSPM_COUPLINGMANAGER_STAGGERED_TPFA_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMASSPM_COUPLINGMANAGER_STAGGERED_TPFA_HH

#include <utility>
#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/boundary/darcydarcy/couplingmapper.hh>

namespace Dumux {

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Coupling manager for Stokes and Darcy domains with equal dimension.
 *        Specialization for staggered-cctpfa coupling.
 */
template<class MDTraits>
class FFMassPMCouplingManagerStaggeredCCTpfa
: public CouplingManager<MDTraits>
{
    using Scalar = typename MDTraits::Scalar;
    using ParentType = CouplingManager<MDTraits>;

public:
    static constexpr auto freeFlowMassIndex = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto porousMediumIndex = typename MDTraits::template SubDomain<1>::Index();

    using SolutionVectorStorage = typename ParentType::SolutionVectorStorage;
private:

    // obtain the type tags of the sub problems
    using FreeFlowMassTypeTag = typename MDTraits::template SubDomain<freeFlowMassIndex>::TypeTag;
    using PorousMediumTypeTag = typename MDTraits::template SubDomain<porousMediumIndex>::TypeTag;

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
        Element<porousMediumIndex> element;
        VolumeVariables<porousMediumIndex> volVars;
        FVElementGeometry<porousMediumIndex> fvGeometry;
        std::size_t dofIdx;
        std::size_t freeFlowMassScvfIdx;
        std::size_t porousMediumScvfIdx;
        mutable VelocityVector velocity; // velocity needs to be set externally, not available in this class
    };

    struct PorousMediumCouplingContext
    {
        Element<freeFlowMassIndex> element;
        VolumeVariables<freeFlowMassIndex> volVars;
        FVElementGeometry<freeFlowMassIndex> fvGeometry;
        std::size_t dofIdx;
        std::size_t porousMediumScvfIdx;
        std::size_t freeFlowMassScvfIdx;
        mutable VelocityVector velocity; // velocity needs to be set externally, not available in this class
    };

    using CouplingMapper = DarcyDarcyBoundaryCouplingMapper<MDTraits>; // TODO rename/generalize class

public:

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    //! Initialize the coupling manager
    void init(std::shared_ptr<Problem<freeFlowMassIndex>> freeFlowMassProblem,
              std::shared_ptr<Problem<porousMediumIndex>> darcyProblem,
              SolutionVectorStorage& curSol)
    {
        this->setSubProblems(std::make_tuple(freeFlowMassProblem, darcyProblem));
        this->attachSolution(curSol);

        couplingMapper_.update(*this);
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

        // we need to update all solution-depenent components of the coupling context
        // the dof of domain J has been deflected
        // if domainJ == freeFlowMassIndex: update volvars in the PorousMediumCouplingContext
        // if domainJ == porousMediumIndex: update volvars in the FreeFlowMassCouplingContext
        // as the update is symmetric we only need to write this once
        auto& context = std::get<1-j>(couplingContext_);
        for (auto& c : context)
        {
            if (c.dofIdx == dofIdxGlobalJ)
            {
                const auto elemSol = elementSolution(c.element, this->curSol(domainJ), this->problem(domainJ).gridGeometry());
                const auto& scv = *scvs(c.fvGeometry).begin();
                c.volVars.update(elemSol, this->problem(domainJ), c.element, scv);
            }
        }
    }

    // \}

    /*!
     * \brief Access the coupling context needed for the Stokes domain
     */
    template<std::size_t i>
    const auto& couplingContext(Dune::index_constant<i> domainI,
                                const FVElementGeometry<i>& fvGeometry,
                                const SubControlVolumeFace<i> scvf) const
    {
        auto& contexts = std::get<i>(couplingContext_);

        if (contexts.empty() || couplingContextBoundForElement_[i] != scvf.insideScvIdx())
            bindCouplingContext_(domainI, fvGeometry);

        for (const auto& context : contexts)
        {
            const auto expectedScvfIdx = domainI == freeFlowMassIndex ? context.freeFlowMassScvfIdx : context.porousMediumScvfIdx;
            if (scvf.index() == expectedScvfIdx)
                return context;
        }

        DUNE_THROW(Dune::InvalidStateException, "No coupling context found at scvf " << scvf.center());
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
        return couplingMapper_.couplingStencil(domainI, eIdx, domainJ);
    }

    // \}

    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI, const SubControlVolumeFace<i>& scvf) const
    { return couplingMapper_.isCoupled(domainI, scvf); }

private:
    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    template<std::size_t i>
    bool isCoupledElement_(Dune::index_constant<i> domainI, std::size_t eIdx) const
    { return couplingMapper_.isCoupledElement(domainI, eIdx); }

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
    template<std::size_t i>
    void bindCouplingContext_(Dune::index_constant<i> domainI, const FVElementGeometry<i>& fvGeometry) const
    {
        auto& context = std::get<domainI>(couplingContext_);
        context.clear();

        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(fvGeometry.element());

        // do nothing if the element is not coupled to the other domain
        if (!isCoupledElement_(domainI, eIdx))
            return;

        couplingContextBoundForElement_[domainI] = eIdx;

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (isCoupled(domainI, scvf))
            {
                const auto otherElementIdx = couplingMapper_.outsideElementIndex(domainI, scvf);
                constexpr auto domainJ = Dune::index_constant<1-domainI>();
                const auto& otherGridGeometry = this->problem(domainJ).gridGeometry();
                const auto& otherElement = otherGridGeometry.element(otherElementIdx);
                auto otherFvGeometry = localView(otherGridGeometry).bindElement(otherElement);

                // there is only one scv for TPFA
                context.push_back({
                    otherElement,
                    volVars_(domainJ, otherElement, *std::begin(scvs(otherFvGeometry))),
                    std::move(otherFvGeometry),
                    otherElementIdx,
                    scvf.index(),
                    couplingMapper_.flipScvfIndex(domainI, scvf),
                    VelocityVector{}
                });
            }
        }
    }

    mutable std::tuple<std::vector<FreeFlowMassCouplingContext>, std::vector<PorousMediumCouplingContext>> couplingContext_;
    mutable std::array<std::size_t, 2> couplingContextBoundForElement_;

    CouplingMapper couplingMapper_;
};

} // end namespace Dumux

#endif
