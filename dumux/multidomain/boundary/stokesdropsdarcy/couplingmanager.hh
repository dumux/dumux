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
 * \ingroup StokesDarcyCoupling
 * \copydoc Dumux::StokesDarcyCouplingManager
 */

#ifndef DUMUX_STOKES_DROPS_DARCY_COUPLINGMANAGER_HH
#define DUMUX_STOKES_DROPS_DARCY_COUPLINGMANAGER_HH

#include <utility>
#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/multidomain/staggeredcouplingmanager.hh>
#include <dumux/discretization/staggered/elementsolution.hh>

#include "couplingdata.hh"
#include "couplingmapper_2in1.hh" // TODO exchange if 3domains used

namespace Dumux {

/*!
 * \ingroup StokesDropsDarcyCoupling
 * \brief Coupling manager for Stokes and Darcy domains with equal dimension
 *           and interface domain with lower dimension.
 */
template<class MDTraits>
class StokesDropsDarcyCouplingManager
: public StaggeredCouplingManagerBase<MDTraits, StokesDropsDarcyCouplingManager<MDTraits>>
{
    using Scalar = typename MDTraits::Scalar;
    using ParentType = StaggeredCouplingManagerBase<MDTraits, StokesDropsDarcyCouplingManager<MDTraits>>;

public:
    static constexpr auto stokesCellCenterIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto stokesFaceIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto cellCenterIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto faceIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto stokesIdx = stokesCellCenterIdx;
    static constexpr auto interfaceIdx = typename MDTraits::template SubDomain<2>::Index();
    static constexpr auto darcyIdx = typename MDTraits::template SubDomain<3>::Index();

private:

    using SolutionVector = typename MDTraits::SolutionVector;

    // obtain the type tags of the sub problems
    using StokesTypeTag = typename MDTraits::template SubDomain<0>::TypeTag;
    using InterfaceTypeTag = typename MDTraits::template SubDomain<2>::TypeTag;
    using DarcyTypeTag = typename MDTraits::template SubDomain<3>::TypeTag;

    using StokesIdType = typename MDTraits::template SubDomain<stokesIdx>::Index;
    using InterfaceIdType = typename MDTraits::template SubDomain<interfaceIdx>::Index;
    using DarcyIdType = typename MDTraits::template SubDomain<darcyIdx>::Index;

    using CouplingStencils = std::unordered_map<std::size_t, std::vector<std::size_t> >;
    using CouplingStencil = CouplingStencils::mapped_type;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    static constexpr bool isCompositional = GetPropType<SubDomainTypeTag<0>, Properties::ModelTraits>::numFluidComponents()> 1;

    template<std::size_t id> using GridView = GetPropType<SubDomainTypeTag<id>, Properties::GridView>;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using NumEqVector = GetPropType<SubDomainTypeTag<id>, Properties::NumEqVector>;
    template<std::size_t id> using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id> using FVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using ElementFluxVariablesCache = typename GetPropType<SubDomainTypeTag<id>, Properties::GridFluxVariablesCache>::LocalView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using PrimaryVariables = typename MDTraits::template SubDomain<id>::PrimaryVariables;
    template<std::size_t id> using SubControlVolume = typename FVGridGeometry<id>::SubControlVolume;
    template<std::size_t id> using SubControlVolumeFace  = typename FVElementGeometry<id>::SubControlVolumeFace;

    using VelocityVector = typename Element<stokesIdx>::Geometry::GlobalCoordinate;
    using GlobalPosition = typename Element<interfaceIdx>::Geometry::GlobalCoordinate;

    using CouplingMapper = StokesDropsDarcyCouplingMapper<MDTraits>;

    // the coupling contexts (information from the respective coupled element)
    struct StokesCouplingContext
    {
        Element<interfaceIdx> element;
        FVElementGeometry<interfaceIdx> fvGeometry;
        VolumeVariables<interfaceIdx> volVars;
        std::size_t stokesScvfIdx;
    };

    struct InterfaceCouplingContext
    {
        Element<stokesIdx> stokesElement;
        FVElementGeometry<stokesIdx> stokesFvGeometry;
        VolumeVariables<stokesIdx> stokesVolVars;
        VelocityVector stokesVelocity;
        Element<darcyIdx> darcyElement;
        FVElementGeometry<darcyIdx> darcyFvGeometry;
        ElementVolumeVariables<darcyIdx> darcyElemVolVars;
        VolumeVariables<darcyIdx> darcyVolVars;
        ElementFluxVariablesCache<darcyIdx> darcyElemFluxVarsCache;
        LocalResidual<darcyIdx> darcyLocalResidual;
        std::size_t darcyElemIdx;
        std::size_t darcyScvfIdx;
        std::size_t interfaceElementIdx;
        VelocityVector darcyElementCenter; // TODO remove, needed for velocity at interface if/pm
    };

    struct DarcyDouplingContext
    {
        Element<interfaceIdx> element;
        FVElementGeometry<interfaceIdx> fvGeometry;
        VolumeVariables<interfaceIdx> volVars;
        std::size_t darcyScvfIdx;
        GlobalPosition elementCenter; // TODO remove, needed for velocity at interface if/pm
    };


public:

    using ParentType::couplingStencil;
    using ParentType::updateCouplingContext;
    using ParentType::evalCouplingResidual;
    using CouplingData = StokesDropsDarcyCouplingData<MDTraits, StokesDropsDarcyCouplingManager<MDTraits>>;

    //! Constructor
    StokesDropsDarcyCouplingManager(std::shared_ptr<const FVGridGeometry<stokesIdx>> stokesFvGridGeometry,
                                    std::shared_ptr<const FVGridGeometry<interfaceIdx>> interfaceFvGridGeometry,
                                    std::shared_ptr<const FVGridGeometry<darcyIdx>> darcyFvGridGeometry) : couplingMapper_(*this)
    { }

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    //! Initialize the coupling manager
    void init(std::shared_ptr<const Problem<stokesIdx>> stokesProblem,
              std::shared_ptr<const Problem<interfaceIdx>> interfaceProblem,
              std::shared_ptr<const Problem<darcyIdx>> darcyProblem,
              const SolutionVector& curSol)
    {
        // TODO compare interface gravity too!
//        if(Dune::FloatCmp::ne(stokesProblem->gravity(), darcyProblem->gravity()))
//              DUNE_THROW(Dune::InvalidStateException, "Both models must use the same gravity vector");

        this->setSubProblems(std::make_tuple(stokesProblem, stokesProblem, interfaceProblem, darcyProblem));
        this->curSol() = curSol;
        couplingData_ = std::make_shared<CouplingData>(*this);
        computeStencils();
    }

    //! Prepare the coupling stencils
    void computeStencils()
    {
        couplingMapper_.computeCouplingMapsAndStencils(stokesCellCenterToInterfaceStencils_,
                                                       stokesFaceToInterfaceStencils_,
                                                       darcyToInterfaceStencils_,
                                                       interfaceToStokesCellCenterStencils_,
                                                       interfaceToStokesFaceStencils_,
                                                       interfaceToDarcyStencils_);

        for(auto&& stencil : stokesCellCenterToInterfaceStencils_)
            removeDuplicates_(stencil.second);
        for(auto&& stencil : stokesFaceToInterfaceStencils_)
            removeDuplicates_(stencil.second);
        for(auto&& stencil : darcyToInterfaceStencils_)
            removeDuplicates_(stencil.second);
        for(auto&& stencil : interfaceToStokesCellCenterStencils_)
            removeDuplicates_(stencil.second);
        for(auto&& stencil : interfaceToStokesFaceStencils_)
            removeDuplicates_(stencil.second);
        for(auto&& stencil : interfaceToDarcyStencils_)
            removeDuplicates_(stencil.second);
    }

    /*!
     * \brief Methods to be accessed by the assembly
     */
    // \{

    template<std::size_t i, class Assembler, std::enable_if_t<(i == stokesCellCenterIdx || i == stokesFaceIdx), int> = 0>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<stokesCellCenterIdx>& element, const Assembler& assembler) const
    { bindCouplingContext(domainI, element); }

    // TODO doc me
    template<std::size_t i, std::enable_if_t<(i == stokesCellCenterIdx || i == stokesFaceIdx), int> = 0>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<stokesCellCenterIdx>& element) const
    {
        stokesCouplingContext_.clear();

        const auto stokesElementIdx = this->problem(stokesIdx).fvGridGeometry().elementMapper().index(element);
        boundStokesElemIdx_ = stokesElementIdx;

        // do nothing if the element is not coupled to the other domain
        if(!couplingMapper_.stokesElementToInterfaceElementMap().count(stokesElementIdx))
            return;

        // prepare the coupling context
        const auto& interfaceIndices = couplingMapper_.stokesElementToInterfaceElementMap().at(stokesElementIdx);
        auto interfaceFvGeometry = localView(this->problem(interfaceIdx).fvGridGeometry());

        // TODO everything needed for coupling?
        for(auto&& indices : interfaceIndices)
        {
            const auto& interfaceElement = this->problem(interfaceIdx).fvGridGeometry().boundingBoxTree().entitySet().entity(indices.eIdx);
            interfaceFvGeometry.bindElement(interfaceElement);
            const auto& scv = (*scvs(interfaceFvGeometry).begin());

            const auto interfaceElemSol = elementSolution(interfaceElement, this->curSol()[interfaceIdx], this->problem(interfaceIdx).fvGridGeometry());
            VolumeVariables<interfaceIdx> interfaceVolVars;
            interfaceVolVars.update(interfaceElemSol, this->problem(interfaceIdx), interfaceElement, scv);

            // add the context
            stokesCouplingContext_.push_back({interfaceElement, interfaceFvGeometry, /*indices.flipScvfIdx, */interfaceVolVars, indices.scvfIdx});
        }
    }

    // TODO doc me
    template<class Assembler>
    void bindCouplingContext(Dune::index_constant<interfaceIdx> domainI, const Element<interfaceIdx>& element, const Assembler& assembler) const
    {
        interfaceCouplingContext_.clear();

        const auto interfaceElementIdx = this->problem(interfaceIdx).fvGridGeometry().elementMapper().index(element);
        boundInterfaceElemIdx_ = interfaceElementIdx;

        // do nothing if the element is not coupled to the other domain
        if(!couplingMapper_.interfaceElementMap().count(interfaceElementIdx))
            return;

        // prepare the coupling context for Stokes
//        const auto& stokesElementIndices = couplingMapper_.interfaceElementMap().at(interfaceElementIdx);
//        assert(stokesElementIndices.size() == 1); // TODO add loop if more than 1 coupled element possible
//        const auto stokesIndex = stokesElementIndices[0];
        const auto coupledIndices = couplingMapper_.interfaceElementMap().at(interfaceElementIdx)[0];

        const auto& stokesElement = this->problem(stokesIdx).fvGridGeometry().boundingBoxTree().entitySet().entity(coupledIndices.stokesEIdx);
        auto stokesFvGeometry = localView(this->problem(stokesIdx).fvGridGeometry());
        stokesFvGeometry.bindElement(stokesElement);

        VelocityVector faceVelocity(0.0);
        for(const auto& scvf : scvfs(stokesFvGeometry))
        {
            if(scvf.index() == coupledIndices.stokesScvfIdx)
                faceVelocity[scvf.directionIndex()] = this->curSol()[stokesFaceIdx][scvf.dofIndex()];
        }

        using PriVarsType = typename VolumeVariables<stokesCellCenterIdx>::PrimaryVariables;
        const auto& cellCenterPriVars = this->curSol()[stokesCellCenterIdx][coupledIndices.stokesEIdx];
        const auto elemSol = makeElementSolutionFromCellCenterPrivars<PriVarsType>(cellCenterPriVars);

        VolumeVariables<stokesIdx> stokesVolVars;
        for(const auto& scv : scvs(stokesFvGeometry))
            stokesVolVars.update(elemSol, this->problem(stokesIdx), stokesElement, scv);

        // prepare the coupling context for Darcy
        const auto& darcyElement = this->problem(darcyIdx).fvGridGeometry().boundingBoxTree().entitySet().entity(coupledIndices.darcyEIdx);
        auto darcyFvGeometry = localView(this->problem(darcyIdx).fvGridGeometry());
        darcyFvGeometry.bindElement(darcyElement);
        const auto& scv = (*scvs(darcyFvGeometry).begin());

        const auto darcyElemSol = elementSolution(darcyElement, this->curSol()[darcyIdx], this->problem(darcyIdx).fvGridGeometry());
        VolumeVariables<darcyIdx> darcyVolVars;
        darcyVolVars.update(darcyElemSol, this->problem(darcyIdx), darcyElement, scv);

        darcyFvGeometry.bind(darcyElement); // TODO necessary?
        auto darcyElemVolVars = Assembler::isImplicit() ? localView(assembler.gridVariables(darcyIdx).curGridVolVars())
                                                        : localView(assembler.gridVariables(darcyIdx).prevGridVolVars());
        auto darcyElemFluxVarsCache = localView(assembler.gridVariables(darcyIdx).gridFluxVarsCache());
        auto darcyLocalResidual = assembler.localResidual(darcyIdx);

        // add the context
        interfaceCouplingContext_.push_back({stokesElement, stokesFvGeometry, stokesVolVars, faceVelocity,
                                             darcyElement, darcyFvGeometry, darcyElemVolVars, darcyVolVars, darcyElemFluxVarsCache, darcyLocalResidual, coupledIndices.darcyEIdx, coupledIndices.darcyScvfIdx,
                                             interfaceElementIdx, scv.center()});

    }

    // TODO doc me
    template<class Assembler>
    void bindCouplingContext(Dune::index_constant<darcyIdx> domainI, const Element<darcyIdx>& element, const Assembler& assembler) const
    { bindCouplingContext(domainI, element); }

    void bindCouplingContext(Dune::index_constant<darcyIdx> domainI, const Element<darcyIdx>& element) const
    {
        darcyCouplingContext_.clear();

        const auto darcyElementIdx = this->problem(darcyIdx).fvGridGeometry().elementMapper().index(element);
        boundDarcyElemIdx_ = darcyElementIdx;

        // do nothing if the element is not coupled to the other domain
        if(!couplingMapper_.darcyElementToInterfaceElementMap().count(darcyElementIdx))
            return;

        // prepare the coupling context
        const auto& interfaceIndices = couplingMapper_.darcyElementToInterfaceElementMap().at(darcyElementIdx);
        auto interfaceFvGeometry = localView(this->problem(interfaceIdx).fvGridGeometry());

        for(auto&& indices : interfaceIndices)
        {
            const auto& interfaceElement = this->problem(interfaceIdx).fvGridGeometry().boundingBoxTree().entitySet().entity(indices.eIdx);
            interfaceFvGeometry.bindElement(interfaceElement);

            const auto interfaceElemSol = elementSolution(interfaceElement, this->curSol()[interfaceIdx], this->problem(interfaceIdx).fvGridGeometry());
            VolumeVariables<interfaceIdx> interfaceVolVars;
            const auto& scv = (*scvs(interfaceFvGeometry).begin());
            interfaceVolVars.update(interfaceElemSol, this->problem(interfaceIdx), interfaceElement, scv);

            // add the context
            darcyCouplingContext_.push_back({interfaceElement, interfaceFvGeometry, interfaceVolVars, indices.scvfIdx, scv.center()});
        }
    }

    /*!
     * \brief Update the coupling context for the interface residual w.r.t. interface DOFs
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<interfaceIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<interfaceIdx> domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<interfaceIdx>& priVarsJ,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
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
     * \brief Update the coupling context for the Stokes residual w.r.t. the interface DOFs
     */
    // TODO copied from stokesdarcy/couplingmanager (Stokes wrt Darcy)
    template<std::size_t i, class LocalAssemblerI, std::enable_if_t<(i == stokesCellCenterIdx || i == stokesFaceIdx), int> = 0>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<interfaceIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<interfaceIdx>& priVars,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ] = priVars;

        for (auto& data : stokesCouplingContext_)
        {
            const auto interfaceElemIdx = this->problem(interfaceIdx).fvGridGeometry().elementMapper().index(data.element);

            if(interfaceElemIdx != dofIdxGlobalJ)
                continue;

            const auto interfaceElemSol = elementSolution(data.element, this->curSol()[interfaceIdx], this->problem(interfaceIdx).fvGridGeometry());

            for(const auto& scv : scvs(data.fvGeometry))
                data.volVars.update(interfaceElemSol, this->problem(interfaceIdx), data.element, scv);
        }
    }

    /*!
     * \brief Update the coupling context for the interface residual w.r.t. the Stokes cell-center DOFs
     */
    // TODO copied from stokesdarcy/couplingmanager (Darcy wrt Stokes cc)
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<interfaceIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<stokesCellCenterIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<stokesCellCenterIdx>& priVars,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ] = priVars;

        for (auto& data : interfaceCouplingContext_)
        {
            const auto stokesElemIdx = this->problem(stokesIdx).fvGridGeometry().elementMapper().index(data.stokesElement);

            if(stokesElemIdx != dofIdxGlobalJ)
                continue;

            using PriVarsType = typename VolumeVariables<stokesCellCenterIdx>::PrimaryVariables;
            const auto elemSol = makeElementSolutionFromCellCenterPrivars<PriVarsType>(priVars);

            for(const auto& scv : scvs(data.stokesFvGeometry))
                data.stokesVolVars.update(elemSol, this->problem(stokesIdx), data.stokesElement, scv);
        }
    }

    /*!
     * \brief Update the coupling context for the interface residual w.r.t. the Stokes face DOFs
     */
    // TODO copied from stokesdarcy/couplingmanager (Darcy wrt Stokes face)
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<interfaceIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<stokesFaceIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<stokesFaceIdx>& priVars,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ] = priVars;

        for (auto& data : interfaceCouplingContext_)
        {
            for(const auto& scvf : scvfs(data.stokesFvGeometry))
            {
                if(scvf.dofIndex() == dofIdxGlobalJ)
                    data.stokesVelocity[scvf.directionIndex()] = priVars;
            }
        }
    }

    /*!
     * \brief Update the coupling context for the Darcy residual w.r.t. the interface DOFs
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(DarcyIdType domainI,
                               const LocalAssemblerI& localAssemblerI,
                               InterfaceIdType domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<interfaceIdx>& priVars,
                               int pvdIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ] = priVars;

        // TODO copied from facetcouplingmanager
//        // Since coupling only occurs via the fluxes, the context does not
//        // have to be updated in explicit time discretization schemes, where
//        // they are strictly evaluated on the old time level
//        if (!LocalAssemblerI::isImplicit())
//            return;

        for (auto& data : darcyCouplingContext_)
        {
            const auto interfaceElemIdx = this->problem(interfaceIdx).fvGridGeometry().elementMapper().index(data.element);

            if(interfaceElemIdx != dofIdxGlobalJ)
                continue;

            const auto interfaceElemSol = elementSolution(data.element, this->curSol()[interfaceIdx], this->problem(interfaceIdx).fvGridGeometry());

            for(const auto& scv : scvs(data.fvGeometry))
            {
                data.volVars.update(interfaceElemSol, this->problem(interfaceIdx), data.element, scv);
                // TODO update necessary? segfault! (extrusion factor)
//                data.darcyElemVolVars[dofIdxGlobalJ].update(darcyElemSol, this->problem(darcyIdx), data.darcyElement, scv);
            }
        }


    }

    /*!
     * \brief Update the coupling context for the interface residual w.r.t. the Darcy DOFs
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(InterfaceIdType domainI,
                               const LocalAssemblerI& localAssemblerI,
                               DarcyIdType domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<darcyIdx>& priVars,
                               int pvdIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ] = priVars;

        // TODO copied from facetcouplingmanager
//        // Since coupling only occurs via the fluxes, the context does not
//        // have to be updated in explicit time discretization schemes, where
//        // they are strictly evaluated on the old time level
//        if (!LocalAssemblerI::isImplicit())
//            return;

        for (auto& data : interfaceCouplingContext_)
        {
            const auto darcyElemIdx = this->problem(darcyIdx).fvGridGeometry().elementMapper().index(data.darcyElement);

            if(darcyElemIdx != dofIdxGlobalJ)
                continue;

            const auto darcyElemSol = elementSolution(data.darcyElement, this->curSol()[darcyIdx], this->problem(darcyIdx).fvGridGeometry());

            for(const auto& scv : scvs(data.darcyFvGeometry))
            {
                data.darcyVolVars.update(darcyElemSol, this->problem(darcyIdx), data.darcyElement, scv);
                // TODO update necessary? segfault! (extrusion factor)
//                data.darcyElemVolVars[dofIdxGlobalJ].update(darcyElemSol, this->problem(darcyIdx), data.darcyElement, scv);

            }
        }
    }


    // \}


    // TODO inherit from stokesdarcycouplingmanager
    /*!
     * \brief Returns whether a given free flow scvf is coupled to the other domain
     */
    bool isCoupledEntity(Dune::index_constant<stokesIdx>, const SubControlVolumeFace<stokesFaceIdx>& scvf) const
    {
        return stokesFaceToInterfaceStencils_.count(scvf.dofIndex());
    }

    // TODO inherit from stokesdarcycouplingmanager
    /*!
     * \brief Returns whether a given free flow scvf is coupled to the other domain
     */
    bool isCoupledEntity(Dune::index_constant<darcyIdx>, const SubControlVolumeFace<darcyIdx>& scvf) const
    {
        return couplingMapper_.isCoupledDarcyScvf(scvf.index());
    }

    /*!
     * \brief The coupling stencils
     */
    // \{

    /*!
     * \brief The Stokes cell center coupling stencil w.r.t. interface DOFs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<stokesCellCenterIdx> domainI,
                                           const Element<stokesIdx>& element,
                                           Dune::index_constant<interfaceIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).fvGridGeometry().elementMapper().index(element);
        if(stokesCellCenterToInterfaceStencils_.count(eIdx))
            return stokesCellCenterToInterfaceStencils_.at(eIdx);
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
     * \brief The interface coupling stencil w.r.t. Stokes cell-centered DOFs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<interfaceIdx> domainI,
                                           const Element<interfaceIdx>& element,
                                           Dune::index_constant<stokesCellCenterIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).fvGridGeometry().elementMapper().index(element);
        if(interfaceToStokesCellCenterStencils_.count(eIdx))
            return interfaceToStokesCellCenterStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief The coupling stencil of a Stokes face w.r.t. interface DOFs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<stokesFaceIdx> domainI,
                                           const SubControlVolumeFace<stokesIdx>& scvf,
                                           Dune::index_constant<interfaceIdx> domainJ) const
    {
        const auto faceDofIdx = scvf.dofIndex();
        if(stokesFaceToInterfaceStencils_.count(faceDofIdx))
            return stokesFaceToInterfaceStencils_.at(faceDofIdx);
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
     * \brief The interface coupling stencil w.r.t. Stokes face DOFs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<interfaceIdx> domainI,
                                           const Element<interfaceIdx>& element,
                                           Dune::index_constant<stokesFaceIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).fvGridGeometry().elementMapper().index(element);
        if (interfaceToStokesFaceStencils_.count(eIdx))
            return interfaceToStokesFaceStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief The Darcy coupling stencil w.r.t. interface DOFs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<darcyIdx> domainI,
                                           const Element<darcyIdx>& element,
                                           Dune::index_constant<interfaceIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).fvGridGeometry().elementMapper().index(element);
        if(darcyToInterfaceStencils_.count(eIdx))
            return darcyToInterfaceStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief The interface coupling stencil w.r.t. Darcy DOFs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<interfaceIdx> domainI,
                                           const Element<interfaceIdx>& element,
                                           Dune::index_constant<darcyIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).fvGridGeometry().elementMapper().index(element);
        if(interfaceToDarcyStencils_.count(eIdx))
            return interfaceToDarcyStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    // \}

    // TODO copied from stokesdarcycouplingmanager, inherit?
    /*!
     * \brief Access the coupling data
     */
    const auto& couplingData() const
    {
        return *couplingData_;
    }

    // TODO copied from stokesdarcycouplingmanager, inherit?
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

    /*! TODO
     * \brief Access the coupling context needed for the interface domain
     */
    const auto& interfaceCouplingContext(const Element<interfaceIdx>& element, const SubControlVolume<interfaceIdx>& scv) const
    {
        // TODO !!! bindCouplingContext(interfaceIdx, element, assembler) --> implement without assembler ?!
//        if (interfaceCouplingContext_.empty())
//            bindCouplingContext(interfaceIdx, element);

        for(const auto& context : interfaceCouplingContext_)
        {
            if(scv.elementIndex() == context.interfaceElementIdx)
                return context;
        }

        // TODO element idx
        DUNE_THROW(Dune::InvalidStateException, "No coupling context found for element ... " /*<< scvf.center()*/);
    }

    // TODO copied from stokesdarcycouplingmanager, inherit?
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
    std::shared_ptr<CouplingData> couplingData_;

    std::unordered_map<std::size_t, std::vector<std::size_t> > stokesCellCenterToInterfaceStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > stokesFaceToInterfaceStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > darcyToInterfaceStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > interfaceToStokesCellCenterStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > interfaceToStokesFaceStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > interfaceToDarcyStencils_;
    std::vector<std::size_t> emptyStencil_;

    ////////////////////////////////////////////////////////////////////////////
    //! The coupling context
    ////////////////////////////////////////////////////////////////////////////
    mutable std::vector<StokesCouplingContext> stokesCouplingContext_;
    mutable std::vector<InterfaceCouplingContext> interfaceCouplingContext_;
    mutable std::vector<DarcyDouplingContext> darcyCouplingContext_;

    mutable std::size_t boundStokesElemIdx_;
    mutable std::size_t boundInterfaceElemIdx_;
    mutable std::size_t boundDarcyElemIdx_;

    CouplingMapper couplingMapper_;
};

} // end namespace Dumux

#endif
