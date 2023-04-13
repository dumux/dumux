// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup StokesDarcyCoupling
 * \copydoc Dumux::StokesDarcyCouplingManager
 */

#ifndef DUMUX_STOKES_DARCY_COUPLINGMANAGER_HH
#define DUMUX_STOKES_DARCY_COUPLINGMANAGER_HH

#include <utility>
#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/multidomain/staggeredcouplingmanager.hh>
#include <dumux/discretization/staggered/elementsolution.hh>

#include "couplingdata.hh"
#include "couplingmapper.hh"

namespace Dumux {

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling manager for Stokes and Darcy domains with equal dimension.
 */
template<class MDTraits>
class StokesDarcyCouplingManager
: public StaggeredCouplingManager<MDTraits>
{
    using Scalar = typename MDTraits::Scalar;
    using ParentType = StaggeredCouplingManager<MDTraits>;

public:
    static constexpr auto stokesFaceIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto stokesCellCenterIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto stokesIdx = stokesFaceIdx;
    static constexpr auto darcyIdx = typename MDTraits::template SubDomain<2>::Index();

private:

    using SolutionVector = typename MDTraits::SolutionVector;

    // obtain the type tags of the sub problems
    using StokesTypeTag = typename MDTraits::template SubDomain<0>::TypeTag;
    using DarcyTypeTag = typename MDTraits::template SubDomain<2>::TypeTag;

    using CouplingStencils = std::unordered_map<std::size_t, std::vector<std::size_t> >;
    using CouplingStencil = CouplingStencils::mapped_type;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    static constexpr bool isCompositional = GetPropType<SubDomainTypeTag<0>, Properties::ModelTraits>::numFluidComponents()> 1;

    template<std::size_t id> using GridView = typename GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>::GridView;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using PrimaryVariables = typename MDTraits::template SubDomain<id>::PrimaryVariables;
    template<std::size_t id> using NumEqVector = Dumux::NumEqVector<PrimaryVariables<id>>;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using GridVolumeVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>;
    template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using ElementBoundaryTypes = GetPropType<SubDomainTypeTag<id>, Properties::ElementBoundaryTypes>;
    template<std::size_t id> using ElementFluxVariablesCache = typename GetPropType<SubDomainTypeTag<id>, Properties::GridFluxVariablesCache>::LocalView;
    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using SubControlVolumeFace  = typename FVElementGeometry<id>::SubControlVolumeFace;

    using CellCenterSolutionVector = GetPropType<StokesTypeTag, Properties::CellCenterSolutionVector>;

    using VelocityVector = typename Element<stokesIdx>::Geometry::GlobalCoordinate;

    using CouplingMapper = StokesDarcyCouplingMapper;

    struct StationaryStokesCouplingContext
    {
        Element<darcyIdx> element;
        FVElementGeometry<darcyIdx> fvGeometry;
        std::size_t darcyScvfIdx;
        std::size_t stokesScvfIdx;
        VolumeVariables<darcyIdx> volVars;
    };

    struct StationaryDarcyCouplingContext
    {
        Element<stokesIdx> element;
        FVElementGeometry<stokesIdx> fvGeometry;
        std::size_t stokesScvfIdx;
        std::size_t darcyScvfIdx;
        VelocityVector velocity;
        VolumeVariables<stokesIdx> volVars;
    };
public:

    using CouplingData = StokesDarcyCouplingData<MDTraits, StokesDarcyCouplingManager<MDTraits>>;

    //! Constructor
    StokesDarcyCouplingManager(std::shared_ptr<const GridGeometry<stokesIdx>> stokesFvGridGeometry,
                               std::shared_ptr<const GridGeometry<darcyIdx>> darcyFvGridGeometry)
    { }

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    //! Initialize the coupling manager
    void init(std::shared_ptr<const Problem<stokesIdx>> stokesProblem,
              std::shared_ptr<const Problem<darcyIdx>> darcyProblem,
              const SolutionVector& curSol)
    {
        if (Dune::FloatCmp::ne(stokesProblem->gravity(), darcyProblem->spatialParams().gravity({})))
            DUNE_THROW(Dune::InvalidStateException, "Both models must use the same gravity vector");

        this->setSubProblems(std::make_tuple(stokesProblem, stokesProblem, darcyProblem));
        this->updateSolution(curSol);
        couplingData_ = std::make_shared<CouplingData>(*this);
        computeStencils();
    }

    // \}

    //! Prepare the coupling stencils
    void computeStencils()
    {
        couplingMapper_.computeCouplingMapsAndStencils(*this,
                                                       darcyToStokesCellCenterCouplingStencils_,
                                                       darcyToStokesFaceCouplingStencils_,
                                                       stokesCellCenterCouplingStencils_,
                                                       stokesFaceCouplingStencils_);

        for(auto&& stencil : darcyToStokesCellCenterCouplingStencils_)
            removeDuplicates_(stencil.second);
        for(auto&& stencil : darcyToStokesFaceCouplingStencils_)
            removeDuplicates_(stencil.second);
        for(auto&& stencil : stokesCellCenterCouplingStencils_)
            removeDuplicates_(stencil.second);
        for(auto&& stencil : stokesFaceCouplingStencils_)
            removeDuplicates_(stencil.second);
    }

    /*!
     * \brief Methods to be accessed by the assembly
     */
    // \{

    using ParentType::evalCouplingResidual;

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Darcy element (i.e. Darcy information)
     */
    template<std::size_t i, class Assembler, std::enable_if_t<(i == stokesCellCenterIdx || i == stokesFaceIdx), int> = 0>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<stokesCellCenterIdx>& element, const Assembler& assembler) const
    { bindCouplingContext(domainI, element); }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Darcy element (i.e. Darcy information)
     */
    template<std::size_t i, std::enable_if_t<(i == stokesCellCenterIdx || i == stokesFaceIdx), int> = 0>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<stokesCellCenterIdx>& element) const
    {
        stokesCouplingContext_.clear();

        const auto stokesElementIdx = this->problem(stokesIdx).gridGeometry().elementMapper().index(element);
        boundStokesElemIdx_ = stokesElementIdx;

        // do nothing if the element is not coupled to the other domain
        if(!couplingMapper_.stokesElementToDarcyElementMap().count(stokesElementIdx))
            return;

        // prepare the coupling context
        const auto& darcyIndices = couplingMapper_.stokesElementToDarcyElementMap().at(stokesElementIdx);
        auto darcyFvGeometry = localView(this->problem(darcyIdx).gridGeometry());
        for(auto&& indices : darcyIndices)
        {
            const auto& darcyElement = this->problem(darcyIdx).gridGeometry().boundingBoxTree().entitySet().entity(indices.eIdx);
            darcyFvGeometry.bindElement(darcyElement);
            const auto& scv = (*scvs(darcyFvGeometry).begin());

            const auto darcyElemSol = elementSolution(darcyElement, this->curSol(darcyIdx), this->problem(darcyIdx).gridGeometry());
            VolumeVariables<darcyIdx> darcyVolVars;
            darcyVolVars.update(darcyElemSol, this->problem(darcyIdx), darcyElement, scv);

            // add the context
            stokesCouplingContext_.push_back({darcyElement, darcyFvGeometry, indices.scvfIdx, indices.flipScvfIdx, darcyVolVars});
        }
    }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Darcy element (i.e. Stokes information)
     */
    template<class Assembler>
    void bindCouplingContext(Dune::index_constant<darcyIdx> domainI, const Element<darcyIdx>& element, const Assembler& assembler) const
    { bindCouplingContext(domainI, element); }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Darcy element (i.e. Stokes information)
     */
    void bindCouplingContext(Dune::index_constant<darcyIdx> domainI, const Element<darcyIdx>& element) const
    {
        darcyCouplingContext_.clear();

        const auto darcyElementIdx = this->problem(darcyIdx).gridGeometry().elementMapper().index(element);
        boundDarcyElemIdx_ = darcyElementIdx;

        // do nothing if the element is not coupled to the other domain
        if(!couplingMapper_.darcyElementToStokesElementMap().count(darcyElementIdx))
            return;

        // prepare the coupling context
        const auto& stokesElementIndices = couplingMapper_.darcyElementToStokesElementMap().at(darcyElementIdx);
        auto stokesFvGeometry = localView(this->problem(stokesIdx).gridGeometry());

        for(auto&& indices : stokesElementIndices)
        {
            const auto& stokesElement = this->problem(stokesIdx).gridGeometry().boundingBoxTree().entitySet().entity(indices.eIdx);
            stokesFvGeometry.bindElement(stokesElement);

            VelocityVector faceVelocity(0.0);

            for(const auto& scvf : scvfs(stokesFvGeometry))
            {
                if(scvf.index() == indices.scvfIdx)
                    faceVelocity[scvf.directionIndex()] = this->curSol(stokesFaceIdx)[scvf.dofIndex()];
            }

            using PriVarsType = typename VolumeVariables<stokesCellCenterIdx>::PrimaryVariables;
            const auto& cellCenterPriVars = this->curSol(stokesCellCenterIdx)[indices.eIdx];
            const auto elemSol = makeElementSolutionFromCellCenterPrivars<PriVarsType>(cellCenterPriVars);

            VolumeVariables<stokesIdx> stokesVolVars;
            for(const auto& scv : scvs(stokesFvGeometry))
                stokesVolVars.update(elemSol, this->problem(stokesIdx), stokesElement, scv);

            // add the context
            darcyCouplingContext_.push_back({stokesElement, stokesFvGeometry, indices.scvfIdx, indices.flipScvfIdx, faceVelocity, stokesVolVars});
        }
    }

    using StaggeredCouplingManager<MDTraits>::updateCouplingContext;

    //! \copydoc StaggeredCouplingManager::updateCouplingContext
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<darcyIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<darcyIdx> domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<darcyIdx>& priVarsJ,
                               int pvIdxJ)
    {
        this->curSol(domainJ)[dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
    }

    //! \copydoc StaggeredCouplingManager::updateCouplingContext
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<darcyIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<stokesCellCenterIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<stokesCellCenterIdx>& priVarsJ,
                               int pvIdxJ)
    {
        this->curSol(domainJ)[dofIdxGlobalJ] = priVarsJ;

        for (auto& data : darcyCouplingContext_)
        {
            const auto stokesElemIdx = this->problem(stokesIdx).gridGeometry().elementMapper().index(data.element);

            if(stokesElemIdx != dofIdxGlobalJ)
                continue;

            using PriVarsType = typename VolumeVariables<stokesCellCenterIdx>::PrimaryVariables;
            const auto elemSol = makeElementSolutionFromCellCenterPrivars<PriVarsType>(priVarsJ);

            for(const auto& scv : scvs(data.fvGeometry))
                data.volVars.update(elemSol, this->problem(stokesIdx), data.element, scv);
        }
    }

    //! \copydoc StaggeredCouplingManager::updateCouplingContext
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<darcyIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<stokesFaceIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<stokesFaceIdx>& priVarsJ,
                               int pvIdxJ)
    {
        this->curSol(domainJ)[dofIdxGlobalJ] = priVarsJ;

        for (auto& data : darcyCouplingContext_)
        {
            for(const auto& scvf : scvfs(data.fvGeometry))
            {
                if(scvf.dofIndex() == dofIdxGlobalJ)
                    data.velocity[scvf.directionIndex()] = priVarsJ;
            }
        }
    }

    //! \copydoc StaggeredCouplingManager::updateCouplingContext
    template<std::size_t i, class LocalAssemblerI, std::enable_if_t<(i == stokesCellCenterIdx || i == stokesFaceIdx), int> = 0>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<darcyIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<darcyIdx>& priVarsJ,
                               int pvIdxJ)
    {
        this->curSol(domainJ)[dofIdxGlobalJ] = priVarsJ;

        for (auto& data : stokesCouplingContext_)
        {
            const auto darcyElemIdx = this->problem(darcyIdx).gridGeometry().elementMapper().index(data.element);

            if(darcyElemIdx != dofIdxGlobalJ)
                continue;

            const auto darcyElemSol = elementSolution(data.element, this->curSol(darcyIdx), this->problem(darcyIdx).gridGeometry());

            for(const auto& scv : scvs(data.fvGeometry))
                data.volVars.update(darcyElemSol, this->problem(darcyIdx), data.element, scv);
        }
    }

    // \}

    /*!
     * \brief Access the coupling data
     */
    const auto& couplingData() const
    {
        return *couplingData_;
    }

    /*!
     * \brief Access the coupling context needed for the Stokes domain
     */
    const auto& stokesCouplingContext(const Element<stokesIdx>& element, const SubControlVolumeFace<stokesIdx>& scvf) const
    {
        if (stokesCouplingContext_.empty() || boundStokesElemIdx_ != scvf.insideScvIdx())
            bindCouplingContext(stokesIdx, element);

        for(const auto& context : stokesCouplingContext_)
        {
            if(scvf.index() == context.stokesScvfIdx)
                return context;
        }

        DUNE_THROW(Dune::InvalidStateException, "No coupling context found at scvf " << scvf.center());
    }

    /*!
     * \brief Access the coupling context needed for the Darcy domain
     */
    const auto& darcyCouplingContext(const Element<darcyIdx>& element, const SubControlVolumeFace<darcyIdx>& scvf) const
    {
        if (darcyCouplingContext_.empty() || boundDarcyElemIdx_ != scvf.insideScvIdx())
            bindCouplingContext(darcyIdx, element);

        for(const auto& context : darcyCouplingContext_)
        {
            if(scvf.index() == context.darcyScvfIdx)
                return context;
        }

        DUNE_THROW(Dune::InvalidStateException, "No coupling context found at scvf " << scvf.center());
    }

    /*!
     * \brief The coupling stencils
     */
    // \{

    using StaggeredCouplingManager<MDTraits>::couplingStencil;

    /*!
     * \brief The Stokes cell center coupling stencil w.r.t. Darcy DOFs.
     *
     * \param domainI The Stokes domain index.
     * \param elementI The Sokes domain element.
     * \param domainJ The Darcy domain index.
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<stokesCellCenterIdx> domainI,
                                           const Element<stokesIdx>& elementI,
                                           Dune::index_constant<darcyIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(elementI);
        if(stokesCellCenterCouplingStencils_.count(eIdx))
            return stokesCellCenterCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief The coupling stencil of domain I, i.e. which domain J DOFs
     *        the given domain I element's residual depends on.
     *
     * \param domainI the index of the domain in which the given element lives.
     * \param elementI the coupled element of domainI
     * \param domainJ the domain index of the coupled domain
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const Element<i>& elementI,
                                           Dune::index_constant<j> domainJ) const
    { return emptyStencil_; }

    /*!
     * \brief Return the Stokes cell indices that influence the residual of
     *        an element in the Darcy domain.
     *
     * \param domainI The darcy domain index.
     * \param elementI The element in the Darcy domain.
     * \param domainJ the domain index of the Stokes domain
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<darcyIdx> domainI,
                                           const Element<darcyIdx>& elementI,
                                           Dune::index_constant<stokesCellCenterIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(elementI);
        if(darcyToStokesCellCenterCouplingStencils_.count(eIdx))
            return darcyToStokesCellCenterCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Return the Stokes face indices that influence the residual of
     *        an element in the Darcy domain.
     *
     * \param domainI The darcy domain index.
     * \param elementI The element in the Darcy domain.
     * \param domainJ the domain index of the Stokes domain
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<darcyIdx> domainI,
                                           const Element<darcyIdx>& elementI,
                                           Dune::index_constant<stokesFaceIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(elementI);
        if (darcyToStokesFaceCouplingStencils_.count(eIdx))
            return darcyToStokesFaceCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Return the dof indices of a subdomain that influence the residual of
     *        a sub-control volume face of the Stokes domain.
     *
     * \param domainI the index of the domain in which the given element lives.
     * \param scvfI the coupled sub-control volume face of the Stokes domain
     * \param domainJ the domain index of the coupled domain
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const SubControlVolumeFace<stokesIdx>& scvfI,
                                           Dune::index_constant<j> domainJ) const
    { return emptyStencil_; }

    /*!
     * \brief The coupling stencil of a Stokes face w.r.t. Darcy DOFs.
     *
     * \param domainI the index of the Stokes domain
     * \param scvfI the coupled subcontrolvolume face of the Stokes domain
     * \param domainJ the index of the Darcy domain
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<stokesFaceIdx> domainI,
                                           const SubControlVolumeFace<stokesIdx>& scvfI,
                                           Dune::index_constant<darcyIdx> domainJ) const
    {
        const auto faceDofIdx = scvfI.dofIndex();
        if(stokesFaceCouplingStencils_.count(faceDofIdx))
            return stokesFaceCouplingStencils_.at(faceDofIdx);
        else
            return emptyStencil_;
    }

    // \}

    /*!
     * \brief There are no additional dof dependencies
     */
    template<class IdType>
    const std::vector<std::size_t>& getAdditionalDofDependencies(IdType id, std::size_t stokesElementIdx) const
    { return emptyStencil_; }

    /*!
     * \brief There are no additional dof dependencies
     */
    template<class IdType>
    const std::vector<std::size_t>& getAdditionalDofDependenciesInverse(IdType id, std::size_t darcyElementIdx) const
    { return emptyStencil_; }

    /*!
     * \brief Returns whether a given free flow scvf is coupled to the other domain
     */
    bool isCoupledEntity(Dune::index_constant<stokesIdx>, const SubControlVolumeFace<stokesFaceIdx>& scvf) const
    {
        return stokesFaceCouplingStencils_.count(scvf.dofIndex());
    }

    /*!
     * \brief Returns whether a given free flow scvf is coupled to the other domain
     */
    bool isCoupledEntity(Dune::index_constant<darcyIdx>, const SubControlVolumeFace<darcyIdx>& scvf) const
    {
        return couplingMapper_.isCoupledDarcyScvf(scvf.index());
    }

protected:

    //! Return a reference to an empty stencil
    std::vector<std::size_t>& emptyStencil()
    { return emptyStencil_; }

    void removeDuplicates_(std::vector<std::size_t>& stencil)
    {
        std::sort(stencil.begin(), stencil.end());
        stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
    }

private:

    std::vector<bool> isCoupledDarcyDof_;
    std::shared_ptr<CouplingData> couplingData_;

    std::unordered_map<std::size_t, std::vector<std::size_t> > stokesCellCenterCouplingStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > stokesFaceCouplingStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > darcyToStokesCellCenterCouplingStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > darcyToStokesFaceCouplingStencils_;
    std::vector<std::size_t> emptyStencil_;

    ////////////////////////////////////////////////////////////////////////////
    //! The coupling context
    ////////////////////////////////////////////////////////////////////////////
    mutable std::vector<StationaryStokesCouplingContext> stokesCouplingContext_;
    mutable std::vector<StationaryDarcyCouplingContext> darcyCouplingContext_;

    mutable std::size_t boundStokesElemIdx_;
    mutable std::size_t boundDarcyElemIdx_;

    CouplingMapper couplingMapper_;
};

} // end namespace Dumux

#endif
