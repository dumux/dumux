// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Freeflow coupling managers (Navier-Stokes mass-momentum coupling)
 */
#ifndef DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_STAGGERED_HH
#define DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_STAGGERED_HH

#include <memory>
#include <tuple>
#include <vector>
#include <deque>

#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/typetraits.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/assemblytraits.hh>
#include <dumux/discretization/facecentered/staggered/consistentlyorientedgrid.hh>

#include <dumux/parallel/parallel_for.hh>
#include <dumux/assembly/coloring.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief The interface of the coupling manager for free flow systems
 * \note coupling manager the face-centered staggered discretization scheme
 */
template<class Traits>
class FCStaggeredFreeFlowCouplingManager
: public CouplingManager<Traits>
{
    using ParentType = CouplingManager<Traits>;
public:
    static constexpr auto freeFlowMomentumIndex = typename Traits::template SubDomain<0>::Index();
    static constexpr auto freeFlowMassIndex = typename Traits::template SubDomain<1>::Index();

    // this can be used if the coupling manager is used inside a meta-coupling manager (e.g. multi-binary)
    // to manager the solution vector storage outside this class
    using SolutionVectorStorage = typename ParentType::SolutionVectorStorage;
private:
    template<std::size_t id> using SubDomainTypeTag = typename Traits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using ElementSeed = typename GridView<id>::Grid::template Codim<0>::EntitySeed;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolume = typename FVElementGeometry<id>::SubControlVolume;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using GridVariables = typename Traits::template SubDomain<id>::GridVariables;
    template<std::size_t id> using ElementVolumeVariables = typename GridVariables<id>::GridVolumeVariables::LocalView;
    template<std::size_t id> using GridFluxVariablesCache = typename GridVariables<id>::GridFluxVariablesCache;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using VolumeVariables = GetPropType<SubDomainTypeTag<id>, Properties::VolumeVariables>;

    using Scalar = typename Traits::Scalar;
    using SolutionVector = typename Traits::SolutionVector;

    template<std::size_t id>
    using SubSolutionVector
        = std::decay_t<decltype(std::declval<SolutionVector>()[Dune::index_constant<id>()])>;

    template<std::size_t id>
    using ConstSubSolutionVectorPtr = const SubSolutionVector<id>*;

    using PrevSolutionVectorStorage = typename Traits::template Tuple<ConstSubSolutionVectorPtr>;

    using CouplingStencilType = std::vector<std::size_t>;

    using GridVariablesTuple = typename Traits::template TupleOfSharedPtr<GridVariables>;

    using FluidSystem = typename VolumeVariables<freeFlowMassIndex>::FluidSystem;

    using VelocityVector = typename SubControlVolumeFace<freeFlowMassIndex>::GlobalPosition;
    static_assert(std::is_same_v<VelocityVector, typename SubControlVolumeFace<freeFlowMomentumIndex>::GlobalPosition>);

    using MassDiscretizationMethod = typename GridGeometry<freeFlowMassIndex>::DiscretizationMethod;

    // This manager identifies a mass dof with its element: the mass sub-control-volume index,
    // the mass dof index and the element index are used interchangeably throughout (see e.g.
    // cellPressure() and couplingStencil()). That only holds for a cell-centered mass domain.
    static_assert(
        MassDiscretizationMethod{} == DiscretizationMethods::cctpfa
        || MassDiscretizationMethod{} == DiscretizationMethods::ccmpfa,
        "The face-centered staggered free-flow coupling manager requires a cell-centered mass domain"
    );

public:

    static constexpr auto pressureIdx = VolumeVariables<freeFlowMassIndex>::Indices::pressureIdx;

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    //! use as regular coupling manager
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> massProblem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol)
    {
        this->setSubProblems(std::make_tuple(momentumProblem, massProblem));
        gridVariables_ = gridVariables;
        this->updateSolution(curSol);

        computeCouplingStencils_();
    }

    //! use as regular coupling manager in a transient setting
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> massProblem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol,
              const SolutionVector& prevSol)
    {
        init(momentumProblem, massProblem, std::forward<GridVariablesTuple>(gridVariables), curSol);

        Dune::Hybrid::forEach(std::make_index_sequence<Traits::numSubDomains>{}, [&](auto i)
        { std::get<i>(prevSolutions_) = &prevSol[i]; });
    }

    //! use as binary coupling manager in multi model context
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> massProblem,
              GridVariablesTuple&& gridVariables,
              const typename ParentType::SolutionVectorStorage& curSol)
    {
        this->setSubProblems(std::make_tuple(momentumProblem, massProblem));
        gridVariables_ = gridVariables;
        this->attachSolution(curSol);

        computeCouplingStencils_();
    }

    //! use as binary coupling manager in multi model context and for transient problems
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> massProblem,
              GridVariablesTuple&& gridVariables,
              const typename ParentType::SolutionVectorStorage& curSol,
              const PrevSolutionVectorStorage& prevSol)
    {
        init(momentumProblem, massProblem, std::forward<GridVariablesTuple>(gridVariables), curSol);
        prevSolutions_ = prevSol;
    }

    // \}

    using CouplingManager<Traits>::evalCouplingResidual;

    /*!
     * \brief evaluates the element residual of a coupled element of domain i which depends on the variables
     *        at the degree of freedom with index dofIdxGlobalJ of domain j
     *
     * \param domainI the domain index of domain i
     * \param localAssemblerI the local assembler assembling the element residual of an element of domain i
     * \param scvI the sub-control-volume of domain i
     * \param domainJ the domain index of domain j
     * \param dofIdxGlobalJ the index of the degree of freedom of domain j which has an influence on the element residual of domain i
     *
     * \note  the element whose residual is to be evaluated can be retrieved from the local assembler
     *        as localAssemblerI.element() as well as all up-to-date variables and caches.
     * \note  the default implementation evaluates the complete element residual
     *        if only parts (i.e. only certain scvs, or only certain terms of the residual) of the residual are coupled
     *        to dof with index dofIdxGlobalJ the function can be overloaded in the coupling manager
     * \return the element residual
     */
    template<std::size_t j, class LocalAssemblerI>
    decltype(auto) evalCouplingResidual(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        const SubControlVolume<freeFlowMomentumIndex>& scvI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ) const
    {
        const auto& problem = localAssemblerI.problem();
        const auto& element = localAssemblerI.element();
        const auto& fvGeometry = localAssemblerI.fvGeometry();
        const auto& curElemVolVars = localAssemblerI.curElemVolVars();
        const auto& prevElemVolVars = localAssemblerI.prevElemVolVars();
        typename LocalAssemblerI::ElementResidualVector residual(localAssemblerI.element().subEntities(1));
        const auto& localResidual = localAssemblerI.localResidual();

        localResidual.evalSource(residual, problem, element, fvGeometry, curElemVolVars, scvI);

        for (const auto& scvf : scvfs(fvGeometry, scvI))
            localResidual.evalFlux(residual, problem, element, fvGeometry, curElemVolVars, localAssemblerI.elemBcTypes(), localAssemblerI.elemFluxVarsCache(), scvf);

        if (!localAssemblerI.assembler().isStationaryProblem())
        {
            assert(isTransient_());
            localResidual.evalStorage(residual, problem, element, fvGeometry, prevElemVolVars, curElemVolVars, scvI);
        }

        return residual;
    }

    /*!
     * \name member functions concerning the coupling stencils
     */
    // \{

    /*!
     * \brief Returns the pressure at a given _frontal_ sub control volume face.
     */
    Scalar pressure(const Element<freeFlowMomentumIndex>& element,
                    const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                    const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        assert(scvf.isFrontal() && !scvf.isLateral() && !scvf.boundary());
        return this->curSol(freeFlowMassIndex)[fvGeometry.elementIndex()][pressureIdx];
    }

    /*!
     * \brief Returns the pressure at the center of a sub control volume corresponding to a given sub control volume face.
     *        This is used for setting a Dirichlet pressure for the mass model when a fixed pressure for the momentum balance is set at another
     *        boundary. Since the the pressure at the given scvf is solution-dependent and thus unknown a priori, we just use the value
     *        of the interior cell here.
     */
    Scalar cellPressure(const Element<freeFlowMassIndex>& element,
                        const SubControlVolumeFace<freeFlowMassIndex>& scvf) const
    {
        return this->curSol(freeFlowMassIndex)[scvf.insideScvIdx()][pressureIdx];
    }

    /*!
     * \brief Returns the density at a given sub control volume face.
     */
    Scalar density(const Element<freeFlowMomentumIndex>& element,
                   const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                   const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                   const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !isTransient_()));
        const auto& insideMomentumScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideMassVolVars = massVolVars_(insideMomentumScv.elementIndex(), considerPreviousTimeStep);

        if (scvf.boundary())
            return insideMassVolVars.density();

        const auto& outsideMomentumScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& outsideMassVolVars = massVolVars_(outsideMomentumScv.elementIndex(), considerPreviousTimeStep);

        // TODO distance weighting
        return 0.5*(insideMassVolVars.density() + outsideMassVolVars.density());
    }

    auto insideAndOutsideDensity(const Element<freeFlowMomentumIndex>& element,
                                 const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                 const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                                 const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !isTransient_()));
        const auto& insideMomentumScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideMassVolVars = massVolVars_(insideMomentumScv.elementIndex(), considerPreviousTimeStep);

        if (scvf.boundary())
            return std::make_pair(insideMassVolVars.density(), insideMassVolVars.density());

        const auto& outsideMomentumScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& outsideMassVolVars = massVolVars_(outsideMomentumScv.elementIndex(), considerPreviousTimeStep);

        return std::make_pair(insideMassVolVars.density(), outsideMassVolVars.density());
    }

    /*!
     * \brief Returns the density at a given sub control volume.
     */
    Scalar density(const Element<freeFlowMomentumIndex>& element,
                   const SubControlVolume<freeFlowMomentumIndex>& scv,
                   const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !isTransient_()));
        return massVolVars_(scv.elementIndex(), considerPreviousTimeStep).density();
    }

    /*!
     * \brief Returns the pressure at a given sub control volume face.
     */
    Scalar effectiveViscosity(const Element<freeFlowMomentumIndex>& element,
                              const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                              const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        const auto& insideMomentumScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideMassVolVars = massVolVars_(insideMomentumScv.elementIndex(), false);

        if (scvf.boundary())
            return insideMassVolVars.viscosity();

        const auto& outsideMomentumScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& outsideMassVolVars = massVolVars_(outsideMomentumScv.elementIndex(), false);

        // TODO distance weighting
        return 0.5*(insideMassVolVars.viscosity() + outsideMassVolVars.viscosity());
    }

    /*!
     * \brief Returns the pressure at a given sub control volume.
     */
    Scalar effectiveViscosity(const Element<freeFlowMomentumIndex>& element,
                              const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                              const SubControlVolume<freeFlowMomentumIndex>& scv) const
    {
        return massVolVars_(scv.elementIndex(), false).viscosity();
    }

    /*!
     * \brief Returns the velocity at a given sub control volume face.
     */
    VelocityVector faceVelocity(const Element<freeFlowMassIndex>& element,
                                const SubControlVolumeFace<freeFlowMassIndex>& scvf) const
    {
        // TODO: rethink this! Maybe we only need scvJ.dofIndex()
        auto momentumFvGeometry = localView(this->problem(freeFlowMomentumIndex).gridGeometry());
        momentumFvGeometry.bindElement(element);

        // the TPFA scvf index corresponds to the staggered scv index (might need mapping)
        const auto localMomentumScvIdx = massScvfToMomentumScvIdx_(scvf, momentumFvGeometry);
        const auto& scvJ = momentumFvGeometry.scv(localMomentumScvIdx);

        // create a unit normal vector oriented in positive coordinate direction
        typename SubControlVolumeFace<freeFlowMassIndex>::GlobalPosition velocity;
        velocity[scvJ.dofAxis()] = 1.0;

        // create the actual velocity vector
        velocity *= this->curSol(freeFlowMomentumIndex)[scvJ.dofIndex()];

        return velocity;
    }

    /*!
     * \brief The coupling stencil of domain I, i.e. which domain J DOFs
     *        the given domain I scv's residual depends on.
     *
     * \param domainI the domain index of domain i
     * \param elementI the coupled element of domain í
     * \param scvI the sub-control volume of domain i
     * \param domainJ the domain index of domain j
     */
    template<std::size_t j>
    const CouplingStencilType& couplingStencil(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                               const Element<freeFlowMomentumIndex>& elementI,
                                               const SubControlVolume<freeFlowMomentumIndex>& scvI,
                                               Dune::index_constant<j> domainJ) const
    { return emptyStencil_; }

    /*!
     * \brief returns an iterable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     *
     * \param domainI the domain index of domain i
     * \param elementI the coupled element of domain í
     * \param domainJ the domain index of domain j
     *
     * \note  The element residual definition depends on the discretization scheme of domain i
     *        box: a container of the residuals of all sub control volumes
     *        cc : the residual of the (sub) control volume
     *        fem: the residual of the element
     * \note  This function has to be implemented by all coupling managers for all combinations of i and j
     */
    const CouplingStencilType& couplingStencil(Dune::index_constant<freeFlowMassIndex> domainI,
                                              const Element<freeFlowMassIndex>& elementI,
                                              Dune::index_constant<freeFlowMomentumIndex> domainJ) const
    {
        const auto eIdx = this->problem(freeFlowMassIndex).gridGeometry().elementMapper().index(elementI);
        return massAndEnergyToMomentumStencils_[eIdx];
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
    const CouplingStencilType& couplingStencil(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                               const Element<freeFlowMomentumIndex>& elementI,
                                               const SubControlVolume<freeFlowMomentumIndex>& scvI,
                                               Dune::index_constant<freeFlowMassIndex> domainJ) const
    {
        return momentumToMassAndEnergyStencils_[scvI.index()];
    }

    // \}

    /*!
     * \name member functions concerning variable caching for element residual evaluations
     */
    // \{

    //! \copydoc CouplingManager::updateCouplingContext
    template<std::size_t i, std::size_t j, class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<j> domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<j>& priVarsJ,
                               int pvIdxJ)
    {
        this->curSol(domainJ)[dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        if constexpr (domainI == freeFlowMomentumIndex && domainJ == freeFlowMassIndex)
        {
            // Only the globally cached volume variables have to be deflected explicitly. Without
            // caching, massVolVars_() recomputes them from the solution, which was just updated
            // above, so the deflection is picked up automatically.
            if constexpr (ElementVolumeVariables<freeFlowMassIndex>::GridVolumeVariables::cachingEnabled)
            {
                const auto& problem = this->problem(domainJ);
                const auto& gridGeometry = problem.gridGeometry();
                const auto deflectedElement = gridGeometry.element(dofIdxGlobalJ);
                auto elemSol = elementSolution(deflectedElement, this->curSol(domainJ), gridGeometry);

                auto fvGeometry = localView(gridGeometry);
                fvGeometry.bindElement(deflectedElement);
                const auto& scv = *(scvs(fvGeometry).begin());
                assert(scv.dofIndex() == dofIdxGlobalJ);

                // Writing into the grid-wide cache is safe under multithreaded assembly: the
                // deflected mass dofs of an element lie within its own stencil, and the coloring
                // (see Dumux::computeColoring) keeps elements that touch the same neighbors -
                // and hence the same volume variables - in different colors.
                gridVars_(freeFlowMassIndex).curGridVolVars().volVars(scv).update(std::move(elemSol), problem, deflectedElement, scv);
            }
        }
    }

    // \}

    /*!
     * \brief Compute colors for multithreaded assembly
     */
    void computeColorsForAssembly()
    {
        // use coloring of the fc staggered discretization for both domains
        elementSets_ = computeColoring(this->problem(freeFlowMomentumIndex).gridGeometry()).sets;
    }

    /*!
     * \brief Execute assembly kernel in parallel
     *
     * \param domainI the domain index of domain i
     * \param assembleElement kernel function to execute for one element
     */
    template<std::size_t i, class AssembleElementFunc>
    void assembleMultithreaded(Dune::index_constant<i> domainI, AssembleElementFunc&& assembleElement) const
    {
        if (elementSets_.empty())
            DUNE_THROW(Dune::InvalidStateException, "Call computeColorsForAssembly before assembling in parallel!");

        // make this element loop run in parallel
        // for this we have to color the elements so that we don't get
        // race conditions when writing into the global matrix
        // each color can be assembled using multiple threads
        const auto& grid = this->problem(freeFlowMomentumIndex).gridGeometry().gridView().grid();
        for (const auto& elements : elementSets_)
        {
            Dumux::parallelFor(elements.size(), [&](const std::size_t eIdx)
            {
                const auto element = grid.entity(elements[eIdx]);
                assembleElement(element);
            });
        }
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


    void computeCouplingStencils_()
    {
        // TODO higher order
        const auto& momentumGridGeometry = this->problem(freeFlowMomentumIndex).gridGeometry();
        auto momentumFvGeometry = localView(momentumGridGeometry);
        massAndEnergyToMomentumStencils_.clear();
        massAndEnergyToMomentumStencils_.resize(momentumGridGeometry.gridView().size(0));

        momentumToMassAndEnergyStencils_.clear();
        momentumToMassAndEnergyStencils_.resize(momentumGridGeometry.numScv());

        for (const auto& element : elements(momentumGridGeometry.gridView()))
        {
            const auto eIdx = momentumGridGeometry.elementMapper().index(element);
            momentumFvGeometry.bind(element);
            for (const auto& scv : scvs(momentumFvGeometry))
            {
                massAndEnergyToMomentumStencils_[eIdx].push_back(scv.dofIndex());
                momentumToMassAndEnergyStencils_[scv.index()].push_back(eIdx);

                // extend the stencil for fluids with variable viscosity and density,
                if constexpr (FluidSystem::isCompressible(0/*phaseIdx*/))
                // if constexpr (FluidSystem::isCompressible(0/*phaseIdx*/) || !FluidSystem::viscosityIsConstant(0/*phaseIdx*/)) // TODO fix on master
                {
                    for (const auto& scvf : scvfs(momentumFvGeometry, scv))
                    {
                        if (scvf.isLateral() && !scvf.boundary())
                        {
                            const auto& outsideScv = momentumFvGeometry.scv(scvf.outsideScvIdx());
                            momentumToMassAndEnergyStencils_[scv.index()].push_back(outsideScv.elementIndex());
                        }
                    }
                }
            }
        }
    }

    std::size_t massScvfToMomentumScvIdx_(const SubControlVolumeFace<freeFlowMassIndex>& massScvf,
                                          const FVElementGeometry<freeFlowMomentumIndex>& momentumFVGeometry) const
    {
        // the faces of the mass balance and the sub-control volumes of the momentum balance are numbered alike,
        // except on a distributed grid, where the mass balance has no faces at the processor boundary
        if (momentumFVGeometry.gridGeometry().gridView().comm().size() == 1)
        {
            if constexpr (ConsistentlyOrientedGrid<typename GridView<freeFlowMomentumIndex>::Grid>{})
                return massScvf.index();
            else
            {
                static const bool makeConsistentlyOriented = getParam<bool>("Grid.MakeConsistentlyOriented", true);
                if (!makeConsistentlyOriented)
                    return massScvf.index();
            }
        }

        for (const auto& momentumScv : scvs(momentumFVGeometry))
        {
            typename SubControlVolumeFace<freeFlowMassIndex>::GlobalPosition momentumUnitOuterNormal(0.0);
            momentumUnitOuterNormal[momentumScv.dofAxis()] = momentumScv.directionSign();
            if (Dune::FloatCmp::eq<typename GridView<freeFlowMomentumIndex>::ctype>(massScvf.unitOuterNormal()*momentumUnitOuterNormal, 1.0))
                return momentumScv.index();
        }
        DUNE_THROW(Dune::InvalidStateException, "No Momentum SCV found");
    }

    CouplingStencilType emptyStencil_;
    std::vector<CouplingStencilType> momentumToMassAndEnergyStencils_;
    std::vector<CouplingStencilType> massAndEnergyToMomentumStencils_;

    /*!
     * \brief Return the mass domain's volume variables of the element with index eIdx
     *
     * There is deliberately no cached coupling context here: with a grid-wide volume variables
     * cache the variables are simply read from it, and without one they are cheap to recompute
     * from the solution. Keeping no mutable state in the manager is what makes the assembly
     * thread-safe (see CouplingManagerSupportsMultithreadedAssembly below).
     *
     * \note With caching enabled this returns a reference into the grid-wide cache, otherwise
     *       a temporary. Bind the result to a const reference to cover both cases.
     */
    decltype(auto) massVolVars_(const std::size_t eIdx, const bool considerPreviousTimeStep) const
    {
        if constexpr (ElementVolumeVariables<freeFlowMassIndex>::GridVolumeVariables::cachingEnabled)
        {
            // for a cell-centered mass domain the scv index equals the element index
            return considerPreviousTimeStep ? gridVars_(freeFlowMassIndex).prevGridVolVars().volVars(eIdx)
                                            : gridVars_(freeFlowMassIndex).curGridVolVars().volVars(eIdx);
        }
        else
        {
            const auto& problem = this->problem(freeFlowMassIndex);
            const auto& gridGeometry = problem.gridGeometry();
            const auto element = gridGeometry.element(eIdx);

            auto fvGeometry = localView(gridGeometry);
            fvGeometry.bindElement(element);
            const auto& scv = *(scvs(fvGeometry).begin());

            const auto& sol = considerPreviousTimeStep ? prevSol_(freeFlowMassIndex)
                                                       : this->curSol(freeFlowMassIndex);

            VolumeVariables<freeFlowMassIndex> volVars;
            volVars.update(elementSolution(element, sol, gridGeometry), problem, element, scv);
            return volVars;
        }
    }

    //! A tuple of std::shared_ptrs to the grid variables of the sub problems
    GridVariablesTuple gridVariables_;

    bool isTransient_() const
    { return std::get<0>(prevSolutions_) != nullptr; }

    template<std::size_t i>
    const SubSolutionVector<i>& prevSol_(Dune::index_constant<i>) const
    { return *std::get<i>(prevSolutions_); }

    PrevSolutionVectorStorage prevSolutions_;

    std::deque<std::vector<ElementSeed<freeFlowMomentumIndex>>> elementSets_;
};

// Multithreaded assembly is supported: the manager keeps no mutable coupling context, so the
// per-element coupling data is derived from the solution independently on each thread. The only
// shared data written during assembly is the mass domain's volume variables cache, which the
// assembly coloring keeps free of races (see updateCouplingContext() above).
template<class T>
struct CouplingManagerSupportsMultithreadedAssembly<FCStaggeredFreeFlowCouplingManager<T>>
: public std::true_type {};

} // end namespace Dumux

#endif
