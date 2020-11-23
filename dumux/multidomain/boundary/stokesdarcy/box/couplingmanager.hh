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

#ifndef DUMUX_STOKES_DARCY_BOX_COUPLINGMANAGER_HH
#define DUMUX_STOKES_DARCY_BOX_COUPLINGMANAGER_HH

#include <utility>
#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/multidomain/staggeredcouplingmanager.hh>
#include <dumux/discretization/staggered/elementsolution.hh>

#include "../couplingmanager.hh"
#include "../couplingdata.hh"
#include "couplingmapper.hh"

namespace Dumux {

namespace Detail {
    enum class ProjectionMethod
    {
        L2Projection, AreaWeightedDofEvaluation
    };

    template <typename T>
    using ProjectionMethodDetector = decltype(std::declval<T>().projectionMethod());

    template<class T>
    static constexpr bool hasProjectionMethod()
    { return Dune::Std::is_detected<ProjectionMethodDetector, T>::value; }


    template<class T>
    static constexpr ProjectionMethod projectionMethod()
    {
        if constexpr (hasProjectionMethod<T>())
            return T::projectionMethod();
        else
            return ProjectionMethod::L2Projection;
    }

    // Each context object contains the data related to one coupling segment
    template <class MDTraits, typename CouplingSegmentGeometry>
    struct StokesCouplingContext
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

        static constexpr auto porousMediumIdx = FreeFlowPorousMediumCouplingManagerBase<MDTraits>::porousMediumIdx;

    public:
        Element<porousMediumIdx> element;
        FVElementGeometry<porousMediumIdx> fvGeometry;
        std::size_t darcyScvfIdx;
        std::size_t stokesScvfIdx;
        std::unique_ptr< ElementVolumeVariables<porousMediumIdx> > elementVolVars;
        std::unique_ptr< ElementFluxVariablesCache<porousMediumIdx> > elementFluxVarsCache;
        CouplingSegmentGeometry segmentGeometry;

        auto permeability() const
        {
            const auto& darcyScvf = fvGeometry.scvf(darcyScvfIdx);
            return (*elementVolVars)[darcyScvf.insideScvIdx()].permeability();
        }
    };

    template<class MDTraits, typename CouplingSegmentGeometry>
    struct DarcyCouplingContext
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

        static constexpr auto freeFlowIdx = FreeFlowPorousMediumCouplingManagerBase<MDTraits>::freeFlowIdx;

        using VelocityVector = typename Element<freeFlowIdx>::Geometry::GlobalCoordinate;

    public:
        Element<freeFlowIdx> element;
        FVElementGeometry<freeFlowIdx> fvGeometry;
        std::size_t stokesScvfIdx;
        std::size_t darcyScvfIdx;
        VelocityVector velocity;
        VolumeVariables<freeFlowIdx> volVars;
        // Darcy context needs information from stokes context because of projection onto stokes faces
        using Container = StokesCouplingContext<MDTraits, CouplingSegmentGeometry>;
        std::unique_ptr< std::vector<Container> > stokesContext;
        CouplingSegmentGeometry segmentGeometry;
    };
};

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling manager for Stokes and Darcy domains with equal dimension.
 */
template<class MDTraits>
class StokesDarcyCouplingManagerImplementation<MDTraits, DiscretizationMethod::box>
: public virtual StaggeredCouplingManager<MDTraits>, public FreeFlowPorousMediumCouplingManagerBase<MDTraits>
{
    using Scalar = typename MDTraits::Scalar;
    using ParentType = StaggeredCouplingManager<MDTraits>;

public:
    static constexpr auto freeFlowFaceIdx = FreeFlowPorousMediumCouplingManagerBase<MDTraits>::freeFlowFaceIdx;
    static constexpr auto freeFlowCellCenterIdx = FreeFlowPorousMediumCouplingManagerBase<MDTraits>::freeFlowCellCenterIdx;
    static constexpr auto freeFlowIdx = FreeFlowPorousMediumCouplingManagerBase<MDTraits>::freeFlowIdx;
    static constexpr auto porousMediumIdx = FreeFlowPorousMediumCouplingManagerBase<MDTraits>::porousMediumIdx;

    using ProjectionMethod = Detail::ProjectionMethod;
    static constexpr auto projectionMethod = Detail::projectionMethod<MDTraits>();
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

    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using NumEqVector = GetPropType<SubDomainTypeTag<id>, Properties::NumEqVector>;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using GridVolumeVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>;
    template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using ElementBoundaryTypes = GetPropType<SubDomainTypeTag<id>, Properties::ElementBoundaryTypes>;
    template<std::size_t id> using ElementFluxVariablesCache = typename GetPropType<SubDomainTypeTag<id>, Properties::GridFluxVariablesCache>::LocalView;
    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using PrimaryVariables = typename MDTraits::template SubDomain<id>::PrimaryVariables;
    template<std::size_t id> using SubControlVolumeFace  = typename FVElementGeometry<id>::SubControlVolumeFace;

    using CellCenterSolutionVector = GetPropType<StokesTypeTag, Properties::CellCenterSolutionVector>;

    using VelocityVector = typename Element<freeFlowIdx>::Geometry::GlobalCoordinate;

    using CouplingMapper = StokesDarcyCouplingMapperBox<MDTraits>;
    using StokesCouplingContext = Detail::StokesCouplingContext<MDTraits, typename CouplingMapper::CouplingSegment::Geometry>;
    using DarcyCouplingContext = Detail::DarcyCouplingContext<MDTraits, typename CouplingMapper::CouplingSegment::Geometry>;

public:

    using ParentType::couplingStencil;
    using ParentType::updateCouplingContext;
    using CouplingData = StokesDarcyCouplingData<MDTraits, StokesDarcyCouplingManager<MDTraits>>;

    //! Constructor
    StokesDarcyCouplingManagerImplementation(std::shared_ptr<const GridGeometry<freeFlowIdx>> stokesFvGridGeometry,
                                             std::shared_ptr<const GridGeometry<porousMediumIdx>> darcyFvGridGeometry)
    { }

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    //! Initialize the coupling manager
    void init(std::shared_ptr<const Problem<freeFlowIdx>> stokesProblem,
              std::shared_ptr<const Problem<porousMediumIdx>> darcyProblem,
              const SolutionVector& curSol)
    {
        if (Dune::FloatCmp::ne(stokesProblem->gravity(), darcyProblem->spatialParams().gravity({})))
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
        couplingMapper_.computeCouplingMapsAndStencils(*this,
                                                       darcyToStokesCellCenterCouplingStencils_,
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

        calculateExtendedStencils_(porousMediumIdx);
    }

    /*!
     * \brief Methods to be accessed by the assembly
     */
    // \{

    using ParentType::evalCouplingResidual;

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Stokes element (i.e. Darcy information)
     */
    template<std::size_t i, class Assembler, std::enable_if_t<(i == freeFlowCellCenterIdx || i == freeFlowFaceIdx), int> = 0>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<freeFlowCellCenterIdx>& element, const Assembler& assembler) const
    {
        stokesCouplingContext_.clear();

        const auto stokesElementIdx = this->problem(freeFlowIdx).gridGeometry().elementMapper().index(element);
        boundStokesElemIdx_ = stokesElementIdx;

        fillCouplingContext(element, assembler, stokesCouplingContext_);
    }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Darcy element (i.e. Stokes information)
     */
    template<class Assembler>
    void bindCouplingContext(Dune::index_constant<porousMediumIdx> domainI, const Element<porousMediumIdx>& element, const Assembler& assembler) const
    {
        darcyCouplingContext_.clear();

        const auto darcyElementIdx = this->problem(porousMediumIdx).gridGeometry().elementMapper().index(element);
        boundDarcyElemIdx_ = darcyElementIdx;

        // do nothing if the element is not coupled to the other domain
        if(!couplingMapper_.darcyElementToStokesElementMap().count(darcyElementIdx))
            return;

        // prepare the coupling context
        const auto& couplingSegments = couplingMapper_.darcyElementToStokesElementMap().at(darcyElementIdx);
        auto stokesFvGeometry = localView(this->problem(freeFlowIdx).gridGeometry());

        for(const auto& couplingSegment : couplingSegments)
        {
            const auto& stokesElement = this->problem(freeFlowIdx).gridGeometry().boundingBoxTree().entitySet().entity(couplingSegment.eIdx);
            stokesFvGeometry.bindElement(stokesElement);

            VelocityVector faceVelocity(0.0);

            for(const auto& scvf : scvfs(stokesFvGeometry))
            {
                if(scvf.index() == couplingSegment.scvfIdx)
                    faceVelocity[scvf.directionIndex()] = this->curSol()[freeFlowFaceIdx][scvf.dofIndex()];
            }

            using PriVarsType = typename VolumeVariables<freeFlowCellCenterIdx>::PrimaryVariables;
            const auto& cellCenterPriVars = this->curSol()[freeFlowCellCenterIdx][couplingSegment.eIdx];
            const auto elemSol = makeElementSolutionFromCellCenterPrivars<PriVarsType>(cellCenterPriVars);

            VolumeVariables<freeFlowIdx> stokesVolVars;
            for(const auto& scv : scvs(stokesFvGeometry))
                stokesVolVars.update(elemSol, this->problem(freeFlowIdx), stokesElement, scv);

            // add the context
            darcyCouplingContext_.push_back({stokesElement, stokesFvGeometry,
                                             couplingSegment.scvfIdx, couplingSegment.flipScvfIdx,
                                             faceVelocity, stokesVolVars,
                                             std::make_unique< std::vector<StokesCouplingContext> >(),
                                             couplingSegment.geometry});

            fillCouplingContext(stokesElement, assembler, (*darcyCouplingContext_.back().stokesContext));
        }
    }

    /*!
     * \brief Update the coupling context for the Darcy residual w.r.t. Darcy DOFs
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<porousMediumIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<porousMediumIdx> domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<porousMediumIdx>& priVarsJ,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
        for (auto& dataI : darcyCouplingContext_)
        {
            const auto& stokesContext = *(dataI.stokesContext);

            for (const auto& dataJ : stokesContext)
            {
                const auto darcyElemSol = elementSolution(dataJ.element, this->curSol()[porousMediumIdx], this->problem(porousMediumIdx).gridGeometry());

                for (auto&& scv : scvs(dataJ.fvGeometry))
                {
                    if(scv.dofIndex() == dofIdxGlobalJ)
                    {
                        auto& volVars = (*dataJ.elementVolVars)[scv];
                        volVars.update(darcyElemSol, this->problem(porousMediumIdx), dataJ.element, scv);
                    }
                }
            }
        }
    }

    /*!
     * \brief Update the coupling context for the Darcy residual w.r.t. the Stokes cell-center DOFs (DarcyToCC)
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<porousMediumIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<freeFlowCellCenterIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<freeFlowCellCenterIdx>& priVars,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ] = priVars;

        for (auto& data : darcyCouplingContext_)
        {
            const auto stokesElemIdx = this->problem(freeFlowIdx).gridGeometry().elementMapper().index(data.element);

            if(stokesElemIdx != dofIdxGlobalJ)
                continue;

            using PriVarsType = typename VolumeVariables<freeFlowCellCenterIdx>::PrimaryVariables;
            const auto elemSol = makeElementSolutionFromCellCenterPrivars<PriVarsType>(priVars);

            for(const auto& scv : scvs(data.fvGeometry))
                data.volVars.update(elemSol, this->problem(freeFlowIdx), data.element, scv);
        }
    }

    /*!
     * \brief Update the coupling context for the Darcy residual w.r.t. the Stokes face DOFs (DarcyToFace)
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<porousMediumIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<freeFlowFaceIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<freeFlowFaceIdx>& priVars,
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
    template<std::size_t i, class LocalAssemblerI, std::enable_if_t<(i == freeFlowCellCenterIdx || i == freeFlowFaceIdx), int> = 0>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<porousMediumIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<porousMediumIdx>& priVars,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ] = priVars;

        for (auto& data : stokesCouplingContext_)
        {
            //Derivatives are assumed to be always calculated with respect to unknowns associated with its own element
            const auto darcyElemSol = elementSolution(data.element, this->curSol()[porousMediumIdx], this->problem(porousMediumIdx).gridGeometry());

            for (auto&& scv : scvs(data.fvGeometry))
            {
                if(scv.dofIndex() == dofIdxGlobalJ)
                {
                    auto& volVars = (*data.elementVolVars)[scv];
                    volVars.update(darcyElemSol, this->problem(porousMediumIdx), data.element, scv);
                }
            }
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
    const auto& stokesCouplingContext(const Element<freeFlowIdx>& element, const SubControlVolumeFace<freeFlowIdx>& scvf) const
    {
        if (stokesCouplingContext_.empty() || boundStokesElemIdx_ != scvf.insideScvIdx())
            DUNE_THROW(Dune::InvalidStateException, "No coupling context found at scvf " << scvf.center());

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
    const auto& darcyCouplingContext(const Element<porousMediumIdx>& element, const SubControlVolumeFace<porousMediumIdx>& scvf) const
    {
        if (darcyCouplingContext_.empty() || boundDarcyElemIdx_ != this->problem(porousMediumIdx).gridGeometry().elementMapper().index(element))
            DUNE_THROW(Dune::InvalidStateException, "No coupling context found at scvf " << scvf.center());

        for(const auto& context : darcyCouplingContext_)
        {
            if(scvf.index() == context.darcyScvfIdx)
                return context;
        }

        DUNE_THROW(Dune::InvalidStateException, "No coupling context found at scvf " << scvf.center());
    }

    /*!
     * \brief Access the coupling context needed for the Darcy domain
     */
    const auto& darcyCouplingContextVector(const Element<porousMediumIdx>& element, const SubControlVolumeFace<porousMediumIdx>& scvf) const
    {
        if (darcyCouplingContext_.empty() || boundDarcyElemIdx_ != this->problem(porousMediumIdx).gridGeometry().elementMapper().index(element))
            DUNE_THROW(Dune::InvalidStateException, "No coupling context found at scvf " << scvf.center());

        return darcyCouplingContext_;
    }

    /*!
     * \brief Access the coupling context needed for the Stokes domain
     */
    const auto& stokesCouplingContextVector(const Element<freeFlowIdx>& element, const SubControlVolumeFace<freeFlowIdx>& scvf) const
    {
        if (stokesCouplingContext_.empty() || boundStokesElemIdx_ != scvf.insideScvIdx())
            DUNE_THROW(Dune::InvalidStateException, "No coupling context found at scvf " << scvf.center());
        return stokesCouplingContext_;
    }

    /*!
     * \brief The coupling stencils
     */
    // \{

    /*!
     * \brief The Stokes cell center coupling stencil w.r.t. Darcy DOFs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<freeFlowCellCenterIdx> domainI,
                                           const Element<freeFlowIdx>& element,
                                           Dune::index_constant<porousMediumIdx> domainJ) const
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
    const CouplingStencil& couplingStencil(Dune::index_constant<porousMediumIdx> domainI,
                                           const Element<porousMediumIdx>& element,
                                           Dune::index_constant<freeFlowCellCenterIdx> domainJ) const
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
    const CouplingStencil& couplingStencil(Dune::index_constant<porousMediumIdx> domainI,
                                           const Element<porousMediumIdx>& element,
                                           Dune::index_constant<freeFlowFaceIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
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
                                           const SubControlVolumeFace<freeFlowIdx>& scvf,
                                           Dune::index_constant<j> domainJ) const
    { return emptyStencil_; }

    /*!
     * \brief The coupling stencil of a Stokes face w.r.t. Darcy DOFs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<freeFlowFaceIdx> domainI,
                                           const SubControlVolumeFace<freeFlowIdx>& scvf,
                                           Dune::index_constant<porousMediumIdx> domainJ) const
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
    bool isCoupledEntity(Dune::index_constant<freeFlowIdx>, const SubControlVolumeFace<freeFlowFaceIdx>& scvf) const
    {
        return stokesFaceCouplingStencils_.count(scvf.dofIndex());
    }

    /*!
     * \brief Returns whether a given free flow scvf is coupled to the other domain
     */
    bool isCoupledEntity(Dune::index_constant<porousMediumIdx>, const Element<porousMediumIdx>& element, const SubControlVolumeFace<porousMediumIdx>& scvf) const
    {
        return couplingMapper_.isCoupledDarcyScvf(this->problem(porousMediumIdx).gridGeometry().elementMapper().index(element),
                                                  scvf.index());
    }

    //! compute the shape function for a given point and geometry
    template<class FVGG, class Geometry, class GlobalPosition, class ShapeValues>
    void getShapeValues(const FVGG& gridGeometry, const Geometry& geo, const GlobalPosition& globalPos, ShapeValues& shapeValues) const
    {
        const auto ipLocal = geo.local(globalPos);
        const auto& localBasis = gridGeometry.feCache().get(geo.type()).localBasis();
        localBasis.evaluateFunction(ipLocal, shapeValues);
    }

     using ParentType::extendJacobianPattern;
    /*!
     * \brief extend the jacobian pattern of the diagonal block of pm domain
     *        by those entries that are not already in the uncoupled pattern
     * \note Such additional dependencies arise from the coupling due to the
     *       projected interface values that are used for diffusive fluxes
     */
    template<class JacobianPattern>
    void extendJacobianPattern(Dune::index_constant<porousMediumIdx> domainI, JacobianPattern& pattern) const
    {
        // add additional dof dependencies
        const auto& gridGeometry = this->problem(domainI).gridGeometry();
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto& dofs = projectionStencil_(domainI, element);

            for (int i = 0; i < element.subEntities(GridView<domainI>::dimension); ++i)
                for (const auto globalJ : dofs)
                    pattern.add(this->problem(domainI).gridGeometry().vertexMapper().subIndex(element, i, GridView<domainI>::dimension), globalJ);
        }
    }

    using ParentType::evalAdditionalDomainDerivatives;
    /*!
     * \brief evaluate additional derivatives of the element residual of a domain with respect
     *        to dofs in the same domain that are not in the regular stencil (per default this is not the case)
     * \note Such additional dependencies can arise from the coupling, e.g. if a coupling source
     *       term depends on a non-local average of a quantity of the same domain
     * \note This is the same for box and cc
     */
    template<class LocalAssemblerI, class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDomainDerivatives(Dune::index_constant<porousMediumIdx> domainI,
                                         const LocalAssemblerI& localAssemblerI,
                                         const typename LocalAssemblerI::LocalResidual::ElementResidualVector&,
                                         JacobianMatrixDiagBlock& A,
                                         GridVariables& gridVariables)
    {
        // Since coupling only occurs via the fluxes, there are no
        // additional derivatives for explicit time integration schemes
        if (!LocalAssemblerI::isImplicit())
            return;

        const auto& curSol = this->curSol()[domainI];
        constexpr auto numEq = std::decay_t<decltype(curSol[0])>::dimension;
        const auto& elementI = localAssemblerI.element();
        const auto& fvGeometryI = localAssemblerI.fvGeometry();

        const auto extendedStencil = extendedStencil_(domainI, elementI);

        // only do something if we have an extended stencil
        if (extendedStencil.empty())
            return;

        auto calcCouplingFluxes = [&](auto& residual)
        {
            for (auto&& scvf : scvfs(fvGeometryI))
            {
                if(isCoupledEntity(domainI, elementI, scvf))
                    localAssemblerI.localResidual().evalFlux(residual,
                                                             this->problem(domainI),
                                                             elementI,
                                                             fvGeometryI,
                                                             localAssemblerI.curElemVolVars(),
                                                             localAssemblerI.elemBcTypes(),
                                                             localAssemblerI.elemFluxVarsCache(),
                                                             scvf);
            }
        };

        // compute the undeflected residual (source only!)
        // initialize the residual vector for all scvs in this element
        typename LocalAssemblerI::LocalResidual::ElementResidualVector origResidual(fvGeometryI.numScv());
        origResidual = 0.0;

        calcCouplingFluxes(origResidual);

        // compute derivate for all additional dofs in the circle stencil
        for (const auto dofIndex : extendedStencil)
        {
            auto partialDerivs = origResidual;
            const auto origPriVars = curSol[dofIndex];

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                // reset partial derivatives
                partialDerivs = 0.0;

                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the coupling context (solution vector and recompute element residual)
                    auto priVars = origPriVars;
                    priVars[pvIdx] = priVar;
                    typename LocalAssemblerI::LocalResidual::ElementResidualVector residual(fvGeometryI.numScv());
                    residual = 0.0;
                    this->updateCouplingContext(domainI, localAssemblerI, domainI, dofIndex, priVars, pvIdx);
                    calcCouplingFluxes(residual);
                    return residual;
                };

                // derive the residuals numerically
                static const int numDiffMethod = getParam<int>("Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, curSol[dofIndex][pvIdx],
                                                          partialDerivs, origResidual, numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (auto&& scvJ : scvs(localAssemblerI.fvGeometry()))
                {
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        // A[i][col][eqIdx][pvIdx] is the rate of change of
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof
                        // 'col'.
                        A[scvJ.dofIndex()][dofIndex][eqIdx][pvIdx] += partialDerivs[scvJ.indexInElement()][eqIdx];
                    }
                }

                // restore the original coupling context
                this->updateCouplingContext(domainI, localAssemblerI, domainI, dofIndex, origPriVars, pvIdx);
            }
        }
    }

protected:

    //! Return a reference to an empty stencil
    std::vector<std::size_t>& emptyStencil()
    { return emptyStencil_; }

    void removeDuplicates_(std::vector<std::size_t>& stencil) const
    {
        std::sort(stencil.begin(), stencil.end());
        stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
    }

private:
            /*!
     * \brief The coupling stencil of a Stokes face w.r.t. Darcy DOFs
     */
    const CouplingStencil& extendedStencil_(Dune::index_constant<porousMediumIdx> domainI,
                                           const Element<porousMediumIdx>& element) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        if (extendedDarcyToDarcyCouplingStencils_.count(eIdx))
            return extendedDarcyToDarcyCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    void calculateExtendedStencils_(Dune::index_constant<porousMediumIdx> domainI)
    {
        extendedDarcyToDarcyCouplingStencils_.clear();
        // add additional dof dependencies
        const auto& gridGeometry = this->problem(domainI).gridGeometry();
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            auto dofs = projectionStencil_(domainI, element);

            auto fvGeometry = localView(gridGeometry);
            fvGeometry.bind(element);
            for (auto&& scv : scvs(fvGeometry))
                dofs.erase(std::remove(dofs.begin(), dofs.end(), scv.dofIndex()), dofs.end());

            auto darcyEIdx = gridGeometry.elementMapper().index(element);
            extendedDarcyToDarcyCouplingStencils_[darcyEIdx] = dofs;
        }
    }

    std::vector<std::size_t> projectionStencil_(Dune::index_constant<porousMediumIdx> domainI, const Element<porousMediumIdx>& element) const
    {
        std::vector<std::size_t> entries;

        const auto& gridGeometry = this->problem(domainI).gridGeometry();
        const auto eIdx = gridGeometry.elementMapper().index(element);
        auto it = darcyToStokesCellCenterCouplingStencils_.find(eIdx);

        // if element is coupled, take one of the neighbors and add coupling stencil to pattern
        if (it != darcyToStokesCellCenterCouplingStencils_.end())
        {
            const auto& couplingStencil = it->second;

            for (auto globalJ : couplingStencil)
            {
                if(stokesCellCenterCouplingStencils_.count(globalJ))
                {
                    const auto& darcyIndices = stokesCellCenterCouplingStencils_.at(globalJ);
                    entries.insert(std::end(entries), std::begin(darcyIndices), std::end(darcyIndices));
                }
            }
        }
        removeDuplicates_(entries);
        return entries;
    }

    template<class Assembler, class CouplingContext>
    void fillCouplingContext(const Element<freeFlowCellCenterIdx>& element, const Assembler& assembler, std::vector<CouplingContext>& couplingContext) const
    {
        const auto stokesElementIdx = this->problem(freeFlowIdx).gridGeometry().elementMapper().index(element);

        // do nothing if the element is not coupled to the other domain
        if(!couplingMapper_.stokesElementToDarcyElementMap().count(stokesElementIdx))
            return;

        // prepare the coupling context
        const auto& couplingSegments = couplingMapper_.stokesElementToDarcyElementMap().at(stokesElementIdx);
        auto darcyFvGeometry = localView(this->problem(porousMediumIdx).gridGeometry());

        for(const auto& couplingSegment : couplingSegments)
        {
            const auto& darcyElement = this->problem(porousMediumIdx).gridGeometry().boundingBoxTree().entitySet().entity(couplingSegment.eIdx);
            darcyFvGeometry.bind(darcyElement);

            auto darcyElemVolVars = localView(assembler.gridVariables(porousMediumIdx).curGridVolVars());
            auto darcyElemFluxVarsCache = localView(assembler.gridVariables(porousMediumIdx).gridFluxVarsCache());

            darcyElemVolVars.bind(darcyElement, darcyFvGeometry, this->curSol()[porousMediumIdx]);
            darcyElemFluxVarsCache.bind(darcyElement, darcyFvGeometry, darcyElemVolVars);

            // add the context
            couplingContext.push_back({darcyElement, darcyFvGeometry, couplingSegment.scvfIdx, couplingSegment.flipScvfIdx,
                                       std::make_unique< ElementVolumeVariables<porousMediumIdx> >( std::move(darcyElemVolVars)),
                                       std::make_unique< ElementFluxVariablesCache<porousMediumIdx> >( std::move(darcyElemFluxVarsCache)),
                                       couplingSegment.geometry});
        }
    }

    std::vector<bool> isCoupledDarcyDof_;
    std::shared_ptr<CouplingData> couplingData_;

    CouplingStencils stokesCellCenterCouplingStencils_;
    CouplingStencils stokesFaceCouplingStencils_;
    CouplingStencils darcyToStokesCellCenterCouplingStencils_;
    std::unordered_map<std::size_t, std::pair<std::vector<std::size_t>,std::vector<std::size_t>> > darcyToStokesFaceCouplingStencils_;
    CouplingStencils extendedDarcyToDarcyCouplingStencils_;
    std::vector<std::size_t> emptyStencil_;

    ////////////////////////////////////////////////////////////////////////////
    //! The coupling context
    ////////////////////////////////////////////////////////////////////////////
    mutable std::vector<StokesCouplingContext> stokesCouplingContext_;
    mutable std::vector<DarcyCouplingContext> darcyCouplingContext_;

    mutable std::size_t boundStokesElemIdx_;
    mutable std::size_t boundDarcyElemIdx_;

    CouplingMapper couplingMapper_;
};

} // end namespace Dumux

#endif
