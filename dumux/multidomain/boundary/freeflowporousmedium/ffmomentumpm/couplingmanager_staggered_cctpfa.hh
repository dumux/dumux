// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FFMomentumPMCouplingManagerStaggeredCCTpfa
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMMOMENTUMPM_COUPLINGMANAGER_STAGGERED_TPFA_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMMOMENTUMPM_COUPLINGMANAGER_STAGGERED_TPFA_HH

#include <utility>
#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/multidomain/couplingmanager.hh>

#include "couplingmapper_staggered_cctpfa.hh"

namespace Dumux {

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Coupling manager for Stokes and Darcy domains with equal dimension
 *        specialization for staggered-cctpfa coupling.
 */
template<class MDTraits>
class FFMomentumPMCouplingManagerStaggeredCCTpfa
: public CouplingManager<MDTraits>
{
    using Scalar = typename MDTraits::Scalar;
    using ParentType = CouplingManager<MDTraits>;

public:
    static constexpr auto freeFlowMomentumIndex = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto porousMediumIndex = typename MDTraits::template SubDomain<1>::Index();

    using SolutionVectorStorage = typename ParentType::SolutionVectorStorage;
private:
    // obtain the type tags of the sub problems
    using FreeFlowMomentumTypeTag = typename MDTraits::template SubDomain<freeFlowMomentumIndex>::TypeTag;
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

    using VelocityVector = typename Element<freeFlowMomentumIndex>::Geometry::GlobalCoordinate;

    struct FreeFlowMomentumCouplingContext
    {
        Element<porousMediumIndex> element;
        VolumeVariables<porousMediumIndex> volVars;
        FVElementGeometry<porousMediumIndex> fvGeometry;
        std::size_t freeFlowMomentumScvfIdx;
        std::size_t porousMediumScvfIdx;
        std::size_t porousMediumDofIdx;
        VelocityVector gravity;
    };

    struct PorousMediumCouplingContext
    {
        Element<freeFlowMomentumIndex> element;
        FVElementGeometry<freeFlowMomentumIndex> fvGeometry;
        std::size_t porousMediumScvfIdx;
        std::size_t freeFlowMomentumScvfIdx;
        std::size_t freeFlowMomentumDofIdx;
        VelocityVector faceVelocity;
    };

    using CouplingMapper = FFMomentumPMCouplingMapperStaggeredCCTpfa<MDTraits, FFMomentumPMCouplingManagerStaggeredCCTpfa<MDTraits>>;

public:

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    //! Initialize the coupling manager
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> freeFlowMomentumProblem,
              std::shared_ptr<Problem<porousMediumIndex>> porousMediumProblem,
              SolutionVectorStorage& curSol)
    {
        this->setSubProblems(std::make_tuple(freeFlowMomentumProblem, porousMediumProblem));
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

        // update the faceVelocity in the PorousMediumCouplingContext
        if constexpr (domainJ == freeFlowMomentumIndex)
        {
            // we only need to update if we are assembling the porous medium domain
            // since the freeflow domain will not use the velocity from the context
            if constexpr (domainI == porousMediumIndex)
            {
                auto& context = std::get<porousMediumIndex>(couplingContext_);
                for (auto& c : context)
                {
                    if (c.freeFlowMomentumDofIdx == dofIdxGlobalJ)
                    {
                        const auto& scvf = localAssemblerI.fvGeometry().scvf(c.porousMediumScvfIdx);
                        c.faceVelocity = faceVelocity(localAssemblerI.element(), scvf);
                    }
                }
            }
        }

        // update volVars in the FreeFlowMomentumCouplingContext
        else if (domainJ == porousMediumIndex)
        {
            auto& context = std::get<freeFlowMomentumIndex>(couplingContext_);
            for (auto& c : context)
            {
                if (c.porousMediumDofIdx == dofIdxGlobalJ)
                {
                    const auto& ggJ = c.fvGeometry.gridGeometry();
                    const auto& scv = *scvs(c.fvGeometry).begin();
                    const auto elemSol = elementSolution(c.element, this->curSol(domainJ), ggJ);
                    c.volVars.update(elemSol, this->problem(domainJ), c.element, scv);
                }
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
            const auto expectedScvfIdx = domainI == freeFlowMomentumIndex ? context.freeFlowMomentumScvfIdx : context.porousMediumScvfIdx;
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
    const CouplingStencil& couplingStencil(Dune::index_constant<porousMediumIndex> domainI,
                                           const Element<porousMediumIndex>& element,
                                           Dune::index_constant<freeFlowMomentumIndex> domainJ) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        return couplingMapper_.couplingStencil(domainI, eIdx, domainJ);
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
                                           Dune::index_constant<porousMediumIndex> domainJ) const
    {
        return couplingMapper_.couplingStencil(domainI, elementI, scvI, domainJ);
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
    { return  couplingMapper_.isCoupledLateralScvf(domainI, scvf); }

    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI, const SubControlVolumeFace<i>& scvf) const
    {
        if constexpr (i == freeFlowMomentumIndex)
            return couplingMapper_.isCoupled(domainI, scvf) || couplingMapper_.isCoupledLateralScvf(domainI, scvf);
        else
            return couplingMapper_.isCoupled(domainI, scvf);
    }

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param scv the sub control volume
     */
    bool isCoupled(Dune::index_constant<freeFlowMomentumIndex> domainI,
                   const SubControlVolume<freeFlowMomentumIndex>& scv) const
    { return couplingMapper_.isCoupled(domainI, scv); }

    /*!
     * \brief Returns the velocity at a given sub control volume face.
     */
    auto faceVelocity(const Element<porousMediumIndex>& element,
                      const SubControlVolumeFace<porousMediumIndex>& scvf) const
    {
        // create a unit normal vector oriented in positive coordinate direction
        auto velocity = scvf.unitOuterNormal();
        using std::abs;
        std::for_each(velocity.begin(), velocity.end(), [](auto& v){ v = abs(v); });

        // create the actual velocity vector
        velocity *= this->curSol(freeFlowMomentumIndex)[
            couplingMapper_.outsideDofIndex(Dune::index_constant<porousMediumIndex>(), scvf)
        ];

        return velocity;
    }

private:
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
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    template<std::size_t i>
    bool isCoupledElement_(Dune::index_constant<i> domainI, std::size_t eIdx) const
    { return couplingMapper_.isCoupledElement(domainI, eIdx); }

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
            if (couplingMapper_.isCoupled(domainI, scvf))
            {
                if constexpr (domainI == freeFlowMomentumIndex)
                {
                    const auto otherElementIdx = couplingMapper_.outsideElementIndex(domainI, scvf);
                    constexpr auto domainJ = Dune::index_constant<1-i>();
                    const auto& otherGridGeometry = this->problem(domainJ).gridGeometry();
                    const auto& otherElement = otherGridGeometry.element(otherElementIdx);
                    auto otherFvGeometry = localView(otherGridGeometry).bindElement(otherElement);

                    // there is only one scv for TPFA
                    context.push_back({
                        otherElement,
                        volVars_(domainJ, otherElement, *std::begin(scvs(otherFvGeometry))),
                        std::move(otherFvGeometry),
                        scvf.index(),
                        couplingMapper_.flipScvfIndex(domainI, scvf),
                        otherElementIdx,
                        this->problem(domainJ).spatialParams().gravity(scvf.center())
                    });
                }

                else if constexpr (domainI == porousMediumIndex)
                {
                    const auto otherElementIdx = couplingMapper_.outsideElementIndex(domainI, scvf);
                    constexpr auto domainJ = Dune::index_constant<1-i>();
                    const auto& otherGridGeometry = this->problem(domainJ).gridGeometry();
                    const auto& otherElement = otherGridGeometry.element(otherElementIdx);
                    auto otherFvGeometry = localView(otherGridGeometry).bindElement(otherElement);

                    const auto otherScvfIdx = couplingMapper_.flipScvfIndex(domainI, scvf);
                    const auto& otherScvf = otherFvGeometry.scvf(otherScvfIdx);
                    const auto& otherScv = otherFvGeometry.scv(otherScvf.insideScvIdx());

                    context.push_back({
                        otherElement,
                        std::move(otherFvGeometry),
                        scvf.index(),
                        otherScvfIdx,
                        otherScv.dofIndex(),
                        faceVelocity(fvGeometry.element(), scvf)
                    });
                }
            }
        }
    }

    mutable std::tuple<std::vector<FreeFlowMomentumCouplingContext>, std::vector<PorousMediumCouplingContext>> couplingContext_;
    mutable std::array<std::size_t, 2> couplingContextBoundForElement_;

    CouplingMapper couplingMapper_;
};

} // end namespace Dumux

#endif
