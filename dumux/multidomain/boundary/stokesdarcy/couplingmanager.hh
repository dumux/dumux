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

#ifndef DUMUX_STOKES_DARCY_COUPLINGMANAGER_HH
#define DUMUX_STOKES_DARCY_COUPLINGMANAGER_HH

#include <utility>
#include <memory>
#include <optional>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
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
    template<std::size_t id> using NumEqVector = GetPropType<SubDomainTypeTag<id>, Properties::NumEqVector>;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using GridVolumeVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>;
    template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using ElementBoundaryTypes = GetPropType<SubDomainTypeTag<id>, Properties::ElementBoundaryTypes>;
    template<std::size_t id> using ElementFluxVariablesCache = typename GetPropType<SubDomainTypeTag<id>, Properties::GridFluxVariablesCache>::LocalView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using PrimaryVariables = typename MDTraits::template SubDomain<id>::PrimaryVariables;
    template<std::size_t id> using SubControlVolumeFace  = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume  = typename FVElementGeometry<id>::SubControlVolume;

    template<std::size_t id> using GridFaceVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridFaceVariables>;
    template<std::size_t id> using ElementFaceVariables = typename GridFaceVariables<id>::LocalView;
    template<std::size_t id> using FaceVariables = typename ElementFaceVariables<id>::FaceVariables;

    template<std::size_t id> using GridVariables = typename MDTraits::template SubDomain<id>::GridVariables;
    using GridVariablesTuple = typename MDTraits::template TupleOfSharedPtr<GridVariables>;

    using CellCenterSolutionVector = GetPropType<StokesTypeTag, Properties::CellCenterSolutionVector>;
    using VelocityVector = typename Element<stokesIdx>::Geometry::GlobalCoordinate;
    using CouplingMapper = StokesDarcyCouplingMapper;
    using TimeLoop = TimeLoopBase<Scalar>;

    struct FreeFlowCouplingContext
    {
        Element<darcyIdx> element;
        FVElementGeometry<darcyIdx> fvGeometry;
        std::size_t darcyScvfIdx;
        std::size_t stokesScvfIdx;
        VolumeVariables<darcyIdx> volVars;

        const auto& getDarcyScvf() const
        { return fvGeometry.scvf(darcyScvfIdx); }
    };

    struct DarcyCouplingContext
    {
        Element<stokesIdx> element;
        FVElementGeometry<stokesIdx> fvGeometry;
        std::size_t stokesScvfIdx;
        std::size_t darcyScvfIdx;

        // quantities only reuqired for CouplingMode::reconstructPorousMediumPressure
        std::optional<VelocityVector> velocity;
        std::optional<VolumeVariables<stokesIdx>> volVars;

        // quantities only reuqired for CouplingMode::reconstructFreeFlowNormalStress
        std::optional<ElementVolumeVariables<stokesIdx>> curElemVolVars;
        std::optional<ElementFaceVariables<stokesIdx>> curElemFaceVars;
        std::optional<ElementVolumeVariables<stokesIdx>> prevElemVolVars;
        std::optional<ElementFaceVariables<stokesIdx>> prevElemFaceVars;
        std::optional<ElementFluxVariablesCache<stokesIdx>> elemFluxVarsCache;

        const auto& getStokesScvf() const
        { return fvGeometry.scvf(stokesScvfIdx); }
    };

public:

    enum class CouplingMode
    {
        reconstructFreeFlowNormalStress, reconstructPorousMediumPressure
    };

    using ParentType::couplingStencil;
    using ParentType::updateCouplingContext;
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
        couplingMode_ = CouplingMode::reconstructPorousMediumPressure;
        printCouplingMode_();

        if (Dune::FloatCmp::ne(stokesProblem->gravity(), darcyProblem->spatialParams().gravity({})))
            DUNE_THROW(Dune::InvalidStateException, "Both models must use the same gravity vector");

        this->setSubProblems(std::make_tuple(stokesProblem, stokesProblem, darcyProblem));
        this->curSol() = curSol;
        couplingData_ = std::make_shared<CouplingData>(*this, timeLoop_);
        computeStencils();
    }

    //! Initialize the coupling manager
    void init(std::shared_ptr<const Problem<stokesIdx>> stokesProblem,
              std::shared_ptr<const Problem<darcyIdx>> darcyProblem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol,
              const CouplingMode couplingMode = CouplingMode::reconstructPorousMediumPressure)
    {
        couplingMode_ = couplingMode;
        gridVariables_ = gridVariables;
        printCouplingMode_();

        if (Dune::FloatCmp::ne(stokesProblem->gravity(), darcyProblem->spatialParams().gravity({})))
            DUNE_THROW(Dune::InvalidStateException, "Both models must use the same gravity vector");

        this->setSubProblems(std::make_tuple(stokesProblem, stokesProblem, darcyProblem));
        this->curSol() = curSol;
        couplingData_ = std::make_shared<CouplingData>(*this, timeLoop_);
        computeStencils();
    }

    //! Initialize the coupling manager
    void init(std::shared_ptr<const Problem<stokesIdx>> stokesProblem,
              std::shared_ptr<const Problem<darcyIdx>> darcyProblem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol,
              const SolutionVector& prevSol,
              std::shared_ptr<TimeLoop> timeLoop,
              const CouplingMode couplingMode = CouplingMode::reconstructPorousMediumPressure)
    {
        timeLoop_ = timeLoop;
        prevSol_ = prevSol;
        init(stokesProblem, darcyProblem, gridVariables, curSol, couplingMode);
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
        couplingMapper_.computeCouplingMapsAndStencils(*this,
                                                       darcyToStokesCellCenterCouplingStencils_,
                                                       darcyToStokesFaceCouplingStencils_,
                                                       stokesCellCenterCouplingStencils_,
                                                       stokesFaceCouplingStencils_,
                                                       couplingMode());

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
        if (!couplingMapper_.stokesElementToDarcyElementMap().count(stokesElementIdx))
            return;

        // prepare the coupling context
        const auto& darcyIndices = couplingMapper_.stokesElementToDarcyElementMap().at(stokesElementIdx);
        auto darcyFvGeometry = localView(this->problem(darcyIdx).gridGeometry());

        for(auto&& indices : darcyIndices)
        {
            const auto& darcyElement = this->problem(darcyIdx).gridGeometry().boundingBoxTree().entitySet().entity(indices.eIdx);
            darcyFvGeometry.bindElement(darcyElement);
            const auto& scv = (*scvs(darcyFvGeometry).begin());

            const auto darcyElemSol = elementSolution(darcyElement, this->curSol()[darcyIdx], this->problem(darcyIdx).gridGeometry());
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
        if (!couplingMapper_.darcyElementToStokesElementMap().count(darcyElementIdx))
            return;

        // prepare the coupling context
        const auto& stokesElementIndices = couplingMapper_.darcyElementToStokesElementMap().at(darcyElementIdx);
        auto stokesFvGeometry = localView(this->problem(stokesIdx).gridGeometry());

        std::optional<VelocityVector> faceVelocity;
        std::optional<VolumeVariables<stokesIdx>> stokesVolVars;
        std::optional<ElementVolumeVariables<stokesIdx>> curStokesElemVolVars;
        std::optional<ElementFaceVariables<stokesIdx>> curStokesElemFaceVars;
        std::optional<ElementVolumeVariables<stokesIdx>> prevStokesElemVolVars;
        std::optional<ElementFaceVariables<stokesIdx>> prevStokesElemFaceVars;
        std::optional<ElementFluxVariablesCache<stokesIdx>> stokesElemFluxVarsCache;

        if (couplingMode() == CouplingMode::reconstructFreeFlowNormalStress)
        {
            curStokesElemVolVars.emplace(localView(gridVars_(stokesIdx).curGridVolVars()));
            curStokesElemFaceVars.emplace(localView(gridVars_(stokesIdx).curGridFaceVars()));
            prevStokesElemVolVars.emplace(localView(gridVars_(stokesIdx).prevGridVolVars()));
            prevStokesElemFaceVars.emplace(localView(gridVars_(stokesIdx).prevGridFaceVars()));
            stokesElemFluxVarsCache.emplace(localView(gridVars_(stokesIdx).gridFluxVarsCache()));
        }

        for (auto&& indices : stokesElementIndices)
        {
            const auto& stokesElement = this->problem(stokesIdx).gridGeometry().boundingBoxTree().entitySet().entity(indices.eIdx);
            stokesFvGeometry.bind(stokesElement);

            if (couplingMode() == CouplingMode::reconstructPorousMediumPressure)
            {
                using PriVarsType = typename VolumeVariables<stokesCellCenterIdx>::PrimaryVariables;
                const auto& cellCenterPriVars = this->curSol()[stokesCellCenterIdx][indices.eIdx];
                const auto elemSol = makeElementSolutionFromCellCenterPrivars<PriVarsType>(cellCenterPriVars);

                stokesVolVars.reset();
                stokesVolVars.emplace();
                for (const auto& scv : scvs(stokesFvGeometry))
                    stokesVolVars->update(elemSol, this->problem(stokesIdx), stokesElement, scv);

                faceVelocity.reset();
                faceVelocity.emplace(0.0);
                for (const auto& scvf : scvfs(stokesFvGeometry))
                {
                    if (scvf.index() == indices.scvfIdx)
                    (*faceVelocity)[scvf.directionIndex()] = this->curSol()[stokesFaceIdx][scvf.dofIndex()];
                }
            }
            else // CouplingMode::reconstructFreeFlowNormalStress
            {
                curStokesElemVolVars->bind(stokesElement, stokesFvGeometry, this->curSol()[stokesCellCenterIdx]);
                curStokesElemFaceVars->bind(stokesElement, stokesFvGeometry, this->curSol()[stokesFaceIdx]);
                stokesElemFluxVarsCache->bind(stokesElement, stokesFvGeometry, *curStokesElemVolVars);
            }

            // add the context
            darcyCouplingContext_.push_back({stokesElement, stokesFvGeometry, indices.scvfIdx, indices.flipScvfIdx,
                                             faceVelocity, stokesVolVars,
                                             curStokesElemVolVars, curStokesElemFaceVars, prevStokesElemVolVars, prevStokesElemFaceVars,
                                             stokesElemFluxVarsCache});
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
            const auto stokesElemIdx = this->problem(stokesIdx).gridGeometry().elementMapper().index(data.element);

            if (stokesElemIdx != dofIdxGlobalJ)
                continue;

            using PriVarsType = typename VolumeVariables<stokesCellCenterIdx>::PrimaryVariables;
            const auto elemSol = makeElementSolutionFromCellCenterPrivars<PriVarsType>(priVars);

            if (couplingMode() == CouplingMode::reconstructPorousMediumPressure)
            {
                for (const auto& scv : scvs(data.fvGeometry))
                    data.volVars->update(elemSol, this->problem(stokesIdx), data.element, scv);
            }
            else // CouplingMode::reconstructFreeFlowNormalStress
            {
                for (const auto& scv : scvs(data.fvGeometry))
                    getVolVarAccess_(domainJ, gridVars_(stokesIdx).curGridVolVars(), data.curElemVolVars.value(), scv).update(elemSol, this->problem(stokesIdx), data.element, scv);
            }
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
                               const PrimaryVariables<stokesFaceIdx>& priVars,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ] = priVars;
        const auto& stokesProblem = this->problem(stokesIdx);

        for (auto& data : darcyCouplingContext_)
        {
            if (couplingMode() == CouplingMode::reconstructPorousMediumPressure)
                (*data.velocity)[data.getStokesScvf().directionIndex()] = priVars;
            else // CouplingMode::reconstructFreeFlowNormalStress
            {
                using FaceSolution = GetPropType<SubDomainTypeTag<stokesIdx>, Properties::StaggeredFaceSolution>;
                FaceSolution faceSol(data.getStokesScvf(), this->curSol()[domainJ], stokesProblem.gridGeometry());
                auto& faceVars = getFaceVarAccess_(domainJ, gridVars_(stokesIdx).curGridFaceVars(), data.curElemFaceVars.value(), data.getStokesScvf());
                faceVars.update(faceSol, stokesProblem, data.element, data.fvGeometry, data.getStokesScvf());
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
            const auto darcyElemIdx = this->problem(darcyIdx).gridGeometry().elementMapper().index(data.element);

            if (darcyElemIdx != dofIdxGlobalJ)
                continue;

            const auto darcyElemSol = elementSolution(data.element, this->curSol()[darcyIdx], this->problem(darcyIdx).gridGeometry());

            for (const auto& scv : scvs(data.fvGeometry))
                data.volVars.update(darcyElemSol, this->problem(darcyIdx), data.element, scv);
        }
    }

    //! Pull up the base class' default implementation
    using ParentType::updateCoupledVariables;

    /*!
     * \brief Update the porous medium flow domain volume variables and flux variables cache
     *        after the coupling context has been updated. This has to be done because the porous medium
              boundary volVars require the updated Dirichlet pressure variable from the free-flow domain.
     */
    template< class LocalAssembler, class UpdatableFluxVarCache >
    void updateCoupledVariables(Dune::index_constant<darcyIdx> domainI,
                                const LocalAssembler& localAssembler,
                                ElementVolumeVariables<darcyIdx>& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        elemVolVars.bind(localAssembler.element(),
                         localAssembler.fvGeometry(),
                         this->curSol()[darcyIdx]);

        // update the transmissibilities subject to the new permeabilities
        elemFluxVarsCache.update(localAssembler.element(),
                                 localAssembler.fvGeometry(),
                                 elemVolVars);
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

        for (const auto& context : stokesCouplingContext_)
        {
            if (scvf.index() == context.stokesScvfIdx)
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

        for (const auto& context : darcyCouplingContext_)
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
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
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
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
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
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        if (darcyToStokesFaceCouplingStencils_.count(eIdx))
            return darcyToStokesFaceCouplingStencils_.at(eIdx);
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

    /*!
     * \brief Returns whether the pointer to the gridVariables was set sucessfully.
     */
    bool gridVariablesAvailable() const
    {
        bool result = true;
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<std::tuple_size_v<GridVariablesTuple>>(), [&](const auto domainId)
        {
            if (!std::get<domainId>(gridVariables_))
                result = false;
        });
        return result;
    }

    /*!
     * \brief Returns the coupling mode.
     */
    CouplingMode couplingMode() const
    { return couplingMode_; }

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

    template<std::size_t i, std::enable_if_t<!getPropValue<SubDomainTypeTag<i>, Properties::EnableGridVolumeVariablesCache>(), int> = 0>
    VolumeVariables<i>& getVolVarAccess_(Dune::index_constant<i> domainIdx, GridVolumeVariables<i>& gridVolVars, ElementVolumeVariables<i>& elemVolVars, const SubControlVolume<i>& scv)
    { return elemVolVars[scv]; }

    template<std::size_t i, std::enable_if_t<getPropValue<SubDomainTypeTag<i>, Properties::EnableGridVolumeVariablesCache>(), int> = 0>
    VolumeVariables<i>& getVolVarAccess_(Dune::index_constant<i> domainIdx, GridVolumeVariables<i>& gridVolVars, ElementVolumeVariables<i>& elemVolVars, const SubControlVolume<i>& scv)
    { return gridVolVars.volVars(scv); }

    template<std::size_t i, std::enable_if_t<!getPropValue<SubDomainTypeTag<i>, Properties::EnableGridFaceVariablesCache>(), int> = 0>
    FaceVariables<i>& getFaceVarAccess_(Dune::index_constant<i> domainIdx, GridFaceVariables<i>& gridFaceVariables, ElementFaceVariables<i>& elemFaceVars, const SubControlVolumeFace<i>& scvf)
    { return elemFaceVars[scvf]; }

    template<std::size_t i, std::enable_if_t<getPropValue<SubDomainTypeTag<i>, Properties::EnableGridFaceVariablesCache>(), int> = 0>
    FaceVariables<i>& getFaceVarAccess_(Dune::index_constant<i> domainIdx, GridFaceVariables<i>& gridFaceVariables, ElementFaceVariables<i>& elemFaceVars, const SubControlVolumeFace<i>& scvf)
    { return gridFaceVariables.faceVars(scvf.index()); }

    void printCouplingMode_() const
    {
        const std::string mode = couplingMode() == CouplingMode::reconstructPorousMediumPressure ? "porous medium pressure" : "free-flow normal stress";
        std::cout << "Coupling mode: reconstruct " << mode << std::endl;
    }


    std::vector<bool> isCoupledDarcyDof_;
    std::shared_ptr<CouplingData> couplingData_;

    std::unordered_map<std::size_t, std::vector<std::size_t> > stokesCellCenterCouplingStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > stokesFaceCouplingStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > darcyToStokesCellCenterCouplingStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > darcyToStokesFaceCouplingStencils_;
    std::vector<std::size_t> emptyStencil_;

    CouplingMode couplingMode_;

    // only for transient problems
    std::shared_ptr<const TimeLoop> timeLoop_ = nullptr;
    const SolutionVector* prevSol_ = nullptr;

    ////////////////////////////////////////////////////////////////////////////
    //! The coupling context
    ////////////////////////////////////////////////////////////////////////////
    mutable std::vector<FreeFlowCouplingContext> stokesCouplingContext_;
    mutable std::vector<DarcyCouplingContext> darcyCouplingContext_;

    mutable std::size_t boundStokesElemIdx_;
    mutable std::size_t boundDarcyElemIdx_;

    /*!
    * \brief A tuple of std::shared_ptrs to the grid variables of the sub problems
    */
    GridVariablesTuple gridVariables_;

    CouplingMapper couplingMapper_;
};

} // end namespace Dumux

#endif
