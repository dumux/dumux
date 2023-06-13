// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FreeFlowMomentumPorousMediumCouplingManager
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMMOMENTUMPM_COUPLINGMANAGER_CVFE_CVFE_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMMOMENTUMPM_COUPLINGMANAGER_CVFE_CVFE_HH

#include <utility>
#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/multidomain/couplingmanager.hh>

#include "couplingmapper_cvfe_cvfe.hh"

namespace Dumux {

namespace DetailFFMomentumPM {
    template<class VolVarType>
    struct VolVarsContainer : public std::vector<VolVarType>
    {
        using VolumeVariables = VolVarType;
    };

    // Each context object contains the data related to one coupling facet
    template <class MDTraits>
    struct FreeFlowMomentumCouplingContext
    {
    private:
        template<std::size_t id>
        using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

        template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
        template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
        template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
        template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
        template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
        template<std::size_t id> using ElementFluxVariablesCache = typename GetPropType<SubDomainTypeTag<id>, Properties::GridFluxVariablesCache>::LocalView;
        template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

        static constexpr auto porousMediumIndex = typename MDTraits::template SubDomain<1>::Index();

    public:
        Element<porousMediumIndex> element;
        FVElementGeometry<porousMediumIndex> fvGeometry;
        std::size_t freeFlowMomentumScvfIdx;
        std::size_t porousMediumScvfIdx;
        VolVarsContainer<VolumeVariables<porousMediumIndex>> elemVolVars;
        // std::size_t facetIdx;

        auto permeability() const
        {
            const auto& pmScvf = fvGeometry.scvf(porousMediumScvfIdx);
            const auto& insideScv = fvGeometry.scv(pmScvf.insideScvIdx());
            return elemVolVars[insideScv.localDofIndex()].permeability();
        }
    };

    template<class MDTraits>
    struct PorousMediumCouplingContext
    {
    private:
        template<std::size_t id>
        using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

        template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
        template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
        template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
        template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
        template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
        template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

        static constexpr auto freeFlowMomentumIndex = typename MDTraits::template SubDomain<0>::Index();

    public:
        Element<freeFlowMomentumIndex> element;
        FVElementGeometry<freeFlowMomentumIndex> fvGeometry;
        std::size_t porousMediumScvfIdx;
        std::size_t freeFlowMomentumScvfIdx;
        VolVarsContainer<VolumeVariables<freeFlowMomentumIndex>> elemVolVars;
        // std::size_t facetIdx;
    };
};

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Coupling manager for Stokes and Darcy domains with equal dimension.
 */
template<class MDTraits>
class FFMomentumPMCouplingManageCvfe
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
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;

    using CouplingMapper = MapperCoupledCvfeSchemes<MDTraits, FFMomentumPMCouplingManageCvfe<MDTraits>>;

    using FreeFlowMomentumCouplingContext = DetailFFMomentumPM::FreeFlowMomentumCouplingContext<MDTraits>;
    using PorousMediumCouplingContext = DetailFFMomentumPM::PorousMediumCouplingContext<MDTraits>;

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

        // update the volVars in the PorousMediumCouplingContext
        if constexpr (domainJ == freeFlowMomentumIndex)
        {
            // we only need to update if we are assembling the porous medium domain
            // since the freeflow domain will not use the context
            if constexpr (domainI == porousMediumIndex)
            {
                auto& context = std::get<porousMediumIndex>(couplingContext_);
                for (auto& c : context)
                {
                    assert(c.fvGeometry.numScv() == c.elemVolVars.size());
                    for (const auto& scv : scvs(c.fvGeometry))
                    {
                        if(scv.dofIndex() == dofIdxGlobalJ)
                        {
                            const auto& ggJ = c.fvGeometry.gridGeometry();
                            const auto elemSol = elementSolution(c.element, this->curSol(domainJ), ggJ);
                            c.elemVolVars[scv.localDofIndex()].update(elemSol, this->problem(domainJ), c.element, scv);
                        }
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
                assert(c.fvGeometry.numScv() == c.elemVolVars.size());
                for (const auto& scv : scvs(c.fvGeometry))
                {
                    if(scv.dofIndex() == dofIdxGlobalJ)
                    {
                        const auto& ggJ = c.fvGeometry.gridGeometry();
                        const auto elemSol = elementSolution(c.element, this->curSol(domainJ), ggJ);
                        c.elemVolVars[scv.localDofIndex()].update(elemSol, this->problem(domainJ), c.element, scv);
                    }
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
     * \brief Stencil of porous domain with respect to free flow momentum domain and vice versa
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const Element<i>& element,
                                           Dune::index_constant<j> domainJ) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        return couplingMapper_.couplingStencil(domainI, eIdx, domainJ);
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

    // /*!
    //  * \brief Returns whether a given scvf is coupled to the other domain
    //  */
    // template<std::size_t i>
    // bool isCoupled(Dune::index_constant<i> domainI, const SubControlVolumeFace<i>& scvf) const
    // {
    //     if constexpr (i == freeFlowMomentumIndex)
    //         return couplingMapper_.isCoupled(domainI, scvf) || couplingMapper_.isCoupledLateralScvf(domainI, scvf);
    //     else
    //         return couplingMapper_.isCoupled(domainI, scvf);
    // }

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param element the element
     * \param scvf the sub control volume face
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI,
                   const Element<i>& element,
                   const SubControlVolumeFace<i>& scvf) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        return couplingMapper_.isCoupled(domainI, eIdx, scvf.index());
    }

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param scv the sub control volume
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI,
                   const SubControlVolume<i>& scv) const
    { return couplingMapper_.isCoupledDof(domainI, scv.dofIndex()); }

    /*!
     * \brief Returns the velocity at a given sub control volume face.
     */
    auto faceVelocity(const Element<porousMediumIndex>& element,
                      const FVElementGeometry<porousMediumIndex>& fvGeometry,
                      const SubControlVolumeFace<porousMediumIndex>& scvf) const
    {
        bindCouplingContext_(Dune::index_constant<porousMediumIndex>(), fvGeometry);
        const auto& context = couplingContext(element, fvGeometry, scvf);

        const auto& ffFvGeometry = context.fvGeometry;
        const auto& localBasis = ffFvGeometry.feLocalBasis();

        std::vector<ShapeValue> shapeValues;
        const auto ipLocal = context.element.geometry().local(scvf.ipGlobal());
        localBasis.evaluateFunction(ipLocal, shapeValues);

        // interpolate velocity at scvf
        VelocityVector velocity(0.0);
        for (const auto& scv : scvs(ffFvGeometry))
            velocity.axpy(shapeValues[scv.localDofIndex()][0], this->curSol(Dune::index_constant<freeFlowMomentumIndex>())[scv.dofIndex()]);

        return velocity;
    }

    // /*!
    //  * \brief Returns the velocity at a given sub control volume face.
    //  */
    // auto faceVelocity(const Element<freeFlowMomentumIndex>& element,
    //                   const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
    //                   const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    // {
    //     const auto& localBasis = fvGeometry.feLocalBasis();

    //     std::vector<ShapeValue> shapeValues;
    //     const auto ipLocal = element.geometry().local(scvf.ipGlobal());
    //     localBasis.evaluateFunction(ipLocal, shapeValues);

    //     // interpolate velocity at scvf
    //     VelocityVector velocity(0.0);
    //     for (const auto& scv : scvs(fvGeometry))
    //         velocity.axpy(shapeValues[scv.localDofIndex()][0], this->curSol(Dune::index_constant<freeFlowMomentumIndex>())[scv.dofIndex()]);

    //     return velocity;
    // }

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
     * \brief prepares all data and variables that are necessary to evaluate the residual depending on coupling terms
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
            if (couplingMapper_.isCoupled(domainI, eIdx, scvf.index()))
            {
                const auto otherElementIdx = couplingMapper_.outsideElementIndex(domainI, eIdx,  scvf.index());
                constexpr auto domainJ = Dune::index_constant<1-i>();
                const auto& otherGridGeometry = this->problem(domainJ).gridGeometry();
                const auto& otherElement = otherGridGeometry.element(otherElementIdx);
                auto otherFvGeometry = localView(otherGridGeometry).bind(otherElement);

                DetailFFMomentumPM::VolVarsContainer<VolumeVariables<1-i>> volVars;
                for (auto&& scv : scvs(otherFvGeometry))
                    volVars.push_back(std::move(volVars_(domainJ, otherElement, scv)));

                context.push_back({
                    otherElement,
                    std::move(otherFvGeometry),
                    scvf.index(),
                    couplingMapper_.flipScvfIndex(domainI, eIdx, scvf.index()),
                    std::move(volVars)
                });
            }
        }
    }

    mutable std::tuple<std::vector<FreeFlowMomentumCouplingContext>, std::vector<PorousMediumCouplingContext>> couplingContext_;
    mutable std::array<std::size_t, 2> couplingContextBoundForElement_;

    CouplingMapper couplingMapper_;
};

} // end namespace Dumux

#endif
