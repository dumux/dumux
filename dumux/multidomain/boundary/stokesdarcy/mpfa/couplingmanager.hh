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
 * \copydoc Dumux::StokesDarcyCouplingManager
 */

#ifndef DUMUX_STOKES_DARCY_COUPLINGMANAGER_HH
#define DUMUX_STOKES_DARCY_COUPLINGMANAGER_HH

#include <utility>
#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/multidomain/staggeredcouplingmanager.hh>

#include "couplingdata.hh"
#include "couplingmapper.hh"

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief Coupling manager for Stokes and Darcy domains with equal dimension.
 */
template<class MDTraits>
class StokesDarcyCouplingManager
: public StaggeredCouplingManagerBase<MDTraits, StokesDarcyCouplingManager<MDTraits>>
{
    using Scalar = typename MDTraits::Scalar;
    using ParentType = StaggeredCouplingManagerBase<MDTraits, StokesDarcyCouplingManager<MDTraits>>;

public:
    static constexpr auto stokesCellCenterIdx = typename MDTraits::template DomainIdx<0>();
    static constexpr auto stokesFaceIdx = typename MDTraits::template DomainIdx<1>();
    static constexpr auto cellCenterIdx = typename MDTraits::template DomainIdx<0>();
    static constexpr auto faceIdx = typename MDTraits::template DomainIdx<1>();
    static constexpr auto stokesIdx = stokesCellCenterIdx;
    static constexpr auto darcyIdx = typename MDTraits::template DomainIdx<2>();

private:

    using SolutionVector = typename MDTraits::SolutionVector;

    // obtain the type tags of the sub problems
    using StokesTypeTag = typename MDTraits::template SubDomainTypeTag<0>;
    using DarcyTypeTag = typename MDTraits::template SubDomainTypeTag<2>;

    using CouplingStencils = std::unordered_map<std::size_t, std::vector<std::size_t> >;
    using CouplingStencil = CouplingStencils::mapped_type;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;

    static constexpr bool isCompositional = GET_PROP_TYPE(SubDomainTypeTag<0>, ModelTraits)::numComponents() > 1;

    template<std::size_t id> using GridView = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridView);
    template<std::size_t id> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);
    template<std::size_t id> using NumEqVector = typename GET_PROP_TYPE(SubDomainTypeTag<id>, NumEqVector);
    template<std::size_t id> using ElementVolumeVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVolumeVariables)::LocalView;
    template<std::size_t id> using GridVolumeVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVolumeVariables);
    template<std::size_t id> using VolumeVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVolumeVariables)::VolumeVariables;
    template<std::size_t id> using FVGridGeometry = typename GET_PROP_TYPE(SubDomainTypeTag<id>, FVGridGeometry);
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using ElementBoundaryTypes = typename GET_PROP_TYPE(SubDomainTypeTag<id>, ElementBoundaryTypes);
    template<std::size_t id> using ElementFluxVariablesCache = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridFluxVariablesCache)::LocalView;
    template<std::size_t id> using GridVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVariables);
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using PrimaryVariables = typename MDTraits::template PrimaryVariables<id>;
    template<std::size_t id> using SubControlVolumeFace  = typename FVElementGeometry<id>::SubControlVolumeFace;

    using CellCenterSolutionVector = typename GET_PROP_TYPE (StokesTypeTag, CellCenterSolutionVector);

    using VelocityVector = typename Element<stokesIdx>::Geometry::GlobalCoordinate;

    using CouplingMapper = StokesDarcyCouplingMapper<MDTraits>;

    struct StationaryStokesCouplingContext
    {
        Element<darcyIdx> element;
        FVElementGeometry<darcyIdx> fvGeometry;
        std::size_t darcyScvfIdx;
        std::size_t stokesScvfIdx;
        VolumeVariables<darcyIdx> volVars;
        std::unique_ptr< ElementVolumeVariables<darcyIdx> > elementVolVars;
        std::unique_ptr< ElementFluxVariablesCache<darcyIdx> > elementFluxVarsCache;
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

    using ParentType::couplingStencil;
    using ParentType::updateCouplingContext;
    using CouplingData = StokesDarcyCouplingData<MDTraits, StokesDarcyCouplingManager<MDTraits>>;

    //! Constructor
    StokesDarcyCouplingManager(std::shared_ptr<const FVGridGeometry<stokesIdx>> stokesFvGridGeometry,
                               std::shared_ptr<const FVGridGeometry<darcyIdx>> darcyFvGridGeometry) : couplingMapper_(*this)
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
        if(Dune::FloatCmp::ne(stokesProblem->gravity(), darcyProblem->gravity()))
            DUNE_THROW(Dune::InvalidStateException, "Both models must use the same gravity vector");

        this->setSubProblems(std::make_tuple(stokesProblem, stokesProblem, darcyProblem));
        this->curSol() = curSol;
        couplingData_ = std::make_shared<CouplingData>(*this);
        computeStencils();
    }

    //! Update after the grid has changed
    void update()
    { }

    // \}

    //! Update the solution vector before assembly
    void updateSolution(const SolutionVector& curSol)
    { this->curSol() = curSol; }

    //! Prepare the coupling stencils
    void computeStencils()
    {
        couplingMapper_.computeCouplingMapsAndStencils(darcyToStokesCellCenterCouplingStencils_,
                                                       darcyToStokesFaceCouplingStencils_,
                                                       stokesCellCenterCouplingStencils_,
                                                       stokesFaceCouplingStencils_);

        for(auto&& stencil : darcyToStokesCellCenterCouplingStencils_)
            removeDuplicates_(stencil.second);
        for(auto&& stencil : darcyToStokesFaceCouplingStencils_)
        {
            removeDuplicates_(stencil.second.first);
            removeDuplicates_(stencil.second.second);
        }
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
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<stokesCellCenterIdx>& element, const Assembler& assembler)
    {
        stokesCouplingContext_.clear();

        const auto stokesElementIdx = this->problem(stokesIdx).fvGridGeometry().elementMapper().index(element);

        // do nothing if the element is not coupled to the other domain
        if(!couplingMapper_.stokesElementToDarcyElementMap().count(stokesElementIdx))
            return;

        // prepare the coupling context
        const auto& darcyIndices = couplingMapper_.stokesElementToDarcyElementMap().at(stokesElementIdx);
        auto darcyFvGeometry = localView(this->problem(darcyIdx).fvGridGeometry());

        for(auto&& indices : darcyIndices)
        {
            const auto& darcyElement = this->problem(darcyIdx).fvGridGeometry().boundingBoxTree().entitySet().entity(indices.eIdx);
            darcyFvGeometry.bind(darcyElement);
            const auto& scv = (*scvs(darcyFvGeometry).begin());

            auto darcyElemVolVars = localView(assembler.gridVariables(darcyIdx).curGridVolVars());
            auto darcyElemFluxVarsCache = localView(assembler.gridVariables(darcyIdx).gridFluxVarsCache());

            darcyElemVolVars.bind(darcyElement, darcyFvGeometry, this->curSol()[darcyIdx]);

            // To bind the darcy coupling context in needed since the darcyElemFluxVarsCache needs the stokes velocities
            static constexpr auto domainId = typename Dune::index_constant<darcyIdx>();
            bindCouplingContext(domainId, darcyElement);
            darcyElemFluxVarsCache.bind(darcyElement, darcyFvGeometry, darcyElemVolVars);

            const auto darcyElemSol = elementSolution(darcyElement, this->curSol()[darcyIdx], this->problem(darcyIdx).fvGridGeometry());
            VolumeVariables<darcyIdx> darcyVolVars;
            darcyVolVars.update(darcyElemSol, this->problem(darcyIdx), darcyElement, scv);

            // add the context
            stokesCouplingContext_.push_back({darcyElement, darcyFvGeometry, indices.scvfIdx, indices.flipScvfIdx, darcyVolVars,
                                              std::make_unique< ElementVolumeVariables<darcyIdx> >( std::move(darcyElemVolVars)),
                                              std::make_unique< ElementFluxVariablesCache<darcyIdx> >( std::move(darcyElemFluxVarsCache))});
        }
    }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Darcy element (i.e. Darcy information)
     */
//    template<std::size_t i, std::enable_if_t<(i == stokesCellCenterIdx || i == stokesFaceIdx), int> = 0>
//    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<stokesCellCenterIdx>& element)
//    {
//        stokesCouplingContext_.clear();
//
//        const auto stokesElementIdx = this->problem(stokesIdx).fvGridGeometry().elementMapper().index(element);
//
//        // do nothing if the element is not coupled to the other domain
//        if(!couplingMapper_.stokesElementToDarcyElementMap().count(stokesElementIdx))
//            return;
//
//        // prepare the coupling context
//        const auto& darcyIndices = couplingMapper_.stokesElementToDarcyElementMap().at(stokesElementIdx);
//        auto darcyFvGeometry = localView(this->problem(darcyIdx).fvGridGeometry());
//
//        for(auto&& indices : darcyIndices)
//        {
//            const auto& darcyElement = this->problem(darcyIdx).fvGridGeometry().boundingBoxTree().entitySet().entity(indices.eIdx);
//            darcyFvGeometry.bindElement(darcyElement);
//            const auto& scv = (*scvs(darcyFvGeometry).begin());
//
//            auto darcyElemVolVars = localView(assembler.gridVariables(bulkId).curGridVolVars());
//            auto darcyElemFluxVarsCache = localView(assembler.gridVariables(bulkId).gridFluxVarsCache());
//
//            const auto darcyElemSol = elementSolution(darcyElement, this->curSol()[darcyIdx], this->problem(darcyIdx).fvGridGeometry());
//            VolumeVariables<darcyIdx> darcyVolVars;
//            darcyVolVars.update(darcyElemSol, this->problem(darcyIdx), darcyElement, scv);
//
//            // add the context
//            stokesCouplingContext_.push_back({darcyElement, darcyFvGeometry, indices.scvfIdx, indices.flipScvfIdx, darcyVolVars});
//        }
//    }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Darcy element (i.e. Stokes information)
     */
    template<class Assembler>
    void bindCouplingContext(Dune::index_constant<darcyIdx> domainI, const Element<darcyIdx>& element, const Assembler& assembler)
    { bindCouplingContext(domainI, element); }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Darcy element (i.e. Stokes information)
     */
    void bindCouplingContext(Dune::index_constant<darcyIdx> domainI, const Element<darcyIdx>& element)
    {
        darcyCouplingContext_.clear();

        const auto darcyElementIdx = this->problem(darcyIdx).fvGridGeometry().elementMapper().index(element);

        // do nothing if the element is not coupled to the other domain
        if(!couplingMapper_.darcyElementToStokesElementMap().count(darcyElementIdx))
            return;

        // prepare the coupling context
        const auto& stokesFaceIndices = darcyToStokesFaceCouplingStencils_.at(darcyElementIdx);

//        if(stokesElementIndices.size() < 2)
//            DUNE_THROW(Dune::InvalidStateException, "The refinement level of the stokes grid should we one level finer");

        const auto& faceDofIdx = stokesFaceIndices.first;
        const auto& scvfIdx = stokesFaceIndices.second;
        auto stokesFvGeometry = localView(this->problem(stokesIdx).fvGridGeometry());

        for(std::size_t i=0; i<scvfIdx.size(); i++)
        {
            const auto& stokesScvf = this->problem(stokesIdx).fvGridGeometry().scvf(scvfIdx[i]);
            const auto& stokesElement = this->problem(stokesIdx).fvGridGeometry().boundingBoxTree().entitySet().entity(stokesScvf.insideScvIdx());
            stokesFvGeometry.bindElement(stokesElement);

            VelocityVector faceVelocity(0.0);
            faceVelocity[stokesScvf.directionIndex()] = this->curSol()[stokesFaceIdx][stokesScvf.dofIndex()];

            using PriVarsType = typename VolumeVariables<stokesCellCenterIdx>::PrimaryVariables;
            const auto& cellCenterPriVars = this->curSol()[stokesCellCenterIdx][stokesScvf.insideScvIdx()];
            const auto elemSol = makeElementSolutionFromCellCenterPrivars<PriVarsType>(cellCenterPriVars);

            VolumeVariables<stokesIdx> stokesVolVars;
            for(const auto& scv : scvs(stokesFvGeometry))
                stokesVolVars.update(elemSol, this->problem(stokesIdx), stokesElement, scv);

            const auto& darcyElementIndices = couplingMapper_.stokesElementToDarcyElementMap().at(stokesScvf.insideScvIdx());
            if(darcyElementIndices.size() > 1)
                DUNE_THROW(Dune::InvalidStateException, "For each stokes element there should be only one darcy element");
            // add the context
            darcyCouplingContext_.push_back({stokesElement, stokesFvGeometry, scvfIdx[i], darcyElementIndices[0].scvfIdx, faceVelocity, stokesVolVars});
        }
    }

    /*!
     * \brief Update the coupling context for the Darcy residual w.r.t. Darcy DOFs
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<darcyIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<darcyIdx> domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<darcyIdx>& priVarsJ,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
    }

    /*!
     * \brief Update the coupling context for the Darcy residual w.r.t. the Stokes cell-center DOFs (DarcyToCC)
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<darcyIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<stokesCellCenterIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<stokesCellCenterIdx>& priVars,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ] = priVars;

        for (auto& data : darcyCouplingContext_)
        {
            const auto stokesElemIdx = this->problem(stokesIdx).fvGridGeometry().elementMapper().index(data.element);

            if(stokesElemIdx != dofIdxGlobalJ)
                continue;

            using PriVarsType = typename VolumeVariables<stokesCellCenterIdx>::PrimaryVariables;
            const auto elemSol = makeElementSolutionFromCellCenterPrivars<PriVarsType>(priVars);

            for(const auto& scv : scvs(data.fvGeometry))
                data.volVars.update(elemSol, this->problem(stokesIdx), data.element, scv);
        }
    }

    /*!
     * \brief Update the coupling context for the Darcy residual w.r.t. the Stokes face DOFs (DarcyToFace)
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<darcyIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<stokesFaceIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<1>& priVars,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ] = priVars;

        for (auto& data : darcyCouplingContext_)
        {
            for(const auto& scvf : scvfs(data.fvGeometry))
            {
                if(scvf.dofIndex() == dofIdxGlobalJ)
                    data.velocity[scvf.directionIndex()] = priVars;
            }
        }
    }

    /*!
     * \brief Update the coupling context for the Stokes cc residual w.r.t. the Darcy DOFs (FaceToDarcy)
     */
    template<std::size_t i, class LocalAssemblerI, std::enable_if_t<(i == stokesCellCenterIdx || i == stokesFaceIdx), int> = 0>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<darcyIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<darcyIdx>& priVars,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ] = priVars;

        for (auto& data : stokesCouplingContext_)
        {
            const auto darcyElemIdx = this->problem(darcyIdx).fvGridGeometry().elementMapper().index(data.element);

            if(darcyElemIdx != dofIdxGlobalJ)
                continue;

            const auto darcyElemSol = elementSolution(data.element, this->curSol()[darcyIdx], this->problem(darcyIdx).fvGridGeometry());

            for(const auto& scv : scvs(data.fvGeometry))
            {
                (*data.elementVolVars)[dofIdxGlobalJ].update(darcyElemSol, this->problem(darcyIdx), data.element, scv);
                data.volVars.update(darcyElemSol, this->problem(darcyIdx), data.element, scv);
            }
            data.elementFluxVarsCache->update(data.element, data.fvGeometry, *data.elementVolVars);
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
    const auto& stokesCouplingContext(const SubControlVolumeFace<stokesIdx>& scvf) const
    {
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
    const auto& darcyCouplingContext(const SubControlVolumeFace<darcyIdx>& scvf) const
    {
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

    /*!
     * \brief The Stokes cell center coupling stencil w.r.t. Darcy DOFs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<stokesCellCenterIdx> domainI,
                                           const Element<stokesIdx>& element,
                                           Dune::index_constant<darcyIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).fvGridGeometry().elementMapper().index(element);
        if(stokesCellCenterCouplingStencils_.count(eIdx))
            return stokesCellCenterCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief The coupling stencil of domain I, i.e. which domain J DOFs
     *        the given domain I element's residual depends on.
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const Element<i>& element,
                                           Dune::index_constant<j> domainJ) const
    { return emptyStencil_; }

    /*!
     * \brief The coupling stencil of domain I, i.e. which domain J dofs
     *        the given domain I element's residual depends on.
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<darcyIdx> domainI,
                                           const Element<darcyIdx>& element,
                                           Dune::index_constant<stokesCellCenterIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).fvGridGeometry().elementMapper().index(element);
        if(darcyToStokesCellCenterCouplingStencils_.count(eIdx))
            return darcyToStokesCellCenterCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief The coupling stencil of domain I, i.e. which domain J dofs
     *        the given domain I element's residual depends on.
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<darcyIdx> domainI,
                                           const Element<darcyIdx>& element,
                                           Dune::index_constant<stokesFaceIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).fvGridGeometry().elementMapper().index(element);
        if (darcyToStokesFaceCouplingStencils_.count(eIdx))
            return darcyToStokesFaceCouplingStencils_.at(eIdx).first;
        else
            return emptyStencil_;
    }

    /*!
     * \brief The coupling stencil of domain I, i.e. which domain J DOFs
     *        the given domain I element's residual depends on.
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const SubControlVolumeFace<stokesIdx>& scvf,
                                           Dune::index_constant<j> domainJ) const
    { return emptyStencil_; }

    /*!
     * \brief The coupling stencil of a Stokes face w.r.t. Darcy DOFs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<stokesFaceIdx> domainI,
                                           const SubControlVolumeFace<stokesIdx>& scvf,
                                           Dune::index_constant<darcyIdx> domainJ) const
    {
        const auto faceDofIdx = scvf.dofIndex();
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
    std::unordered_map<std::size_t, std::pair<std::vector<std::size_t>,std::vector<std::size_t>> > darcyToStokesFaceCouplingStencils_;
    std::vector<std::size_t> emptyStencil_;

    ////////////////////////////////////////////////////////////////////////////
    //! The coupling context
    ////////////////////////////////////////////////////////////////////////////
    std::vector<StationaryStokesCouplingContext> stokesCouplingContext_;
    std::vector<StationaryDarcyCouplingContext> darcyCouplingContext_;

    CouplingMapper couplingMapper_;
};

} // end namespace Dumux

#endif
