// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Freeflow coupling managers (Navier-Stokes mass-momentum coupling)
 */
#ifndef DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_PQ1BUBBLE_HH
#define DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_PQ1BUBBLE_HH

#warning "This file is deprecated and will be removed after 3.7. Use CVFEFreeFlowCouplingManager."

#include <memory>
#include <tuple>
#include <vector>
#include <deque>

#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/typetraits.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/fvassembler.hh>

#include <dumux/parallel/parallel_for.hh>
#include <dumux/assembly/coloring.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief The interface of the coupling manager for free flow systems
 * \note coupling manager for the control-volume finite-element discretization scheme
 */
template<class Traits>
class PQ1BubbleFreeFlowCouplingManager
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

    using CouplingStencilType = std::vector<std::size_t>;

    using GridVariablesTuple = typename Traits::template TupleOfSharedPtr<GridVariables>;

    using FluidSystem = typename VolumeVariables<freeFlowMassIndex>::FluidSystem;

    using VelocityVector = typename SubControlVolumeFace<freeFlowMassIndex>::GlobalPosition;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;

    static_assert(std::is_same_v<VelocityVector, typename SubControlVolumeFace<freeFlowMomentumIndex>::GlobalPosition>);

    struct MomentumCouplingContext
    {
        FVElementGeometry<freeFlowMassIndex> fvGeometry;
        ElementVolumeVariables<freeFlowMassIndex> curElemVolVars;
        ElementVolumeVariables<freeFlowMassIndex> prevElemVolVars;
        std::size_t eIdx;
    };

    struct MassAndEnergyCouplingContext
    {
        MassAndEnergyCouplingContext(FVElementGeometry<freeFlowMomentumIndex>&& f, const std::size_t i)
        : fvGeometry(std::move(f))
        , eIdx(i)
        {}

        FVElementGeometry<freeFlowMomentumIndex> fvGeometry;
        std::size_t eIdx;
    };

    using MomentumDiscretizationMethod = typename GridGeometry<freeFlowMomentumIndex>::DiscretizationMethod;
    using MassDiscretizationMethod = typename GridGeometry<freeFlowMassIndex>::DiscretizationMethod;

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
        prevSol_ = &prevSol;
        isTransient_ = true;
    }

    //! use as binary coupling manager in multi model context
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> massProblem,
              GridVariablesTuple&& gridVariables,
              typename ParentType::SolutionVectorStorage& curSol)
    {
        this->setSubProblems(std::make_tuple(momentumProblem, massProblem));
        gridVariables_ = gridVariables;
        this->attachSolution(curSol);

        computeCouplingStencils_();
    }

    // \}

    /*!
     * \name member functions concerning the coupling stencils
     */
    // \{

    /*!
     * \brief Returns the pressure at a given sub control volume face
     */
    Scalar pressure(const Element<freeFlowMomentumIndex>& element,
                    const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                    const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                    const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !this->isTransient_));
        const auto& gg = this->problem(freeFlowMassIndex).gridGeometry();
        const auto& sol = considerPreviousTimeStep ? (*prevSol_)[freeFlowMassIndex]
                                                   :  this->curSol(freeFlowMassIndex);
        const auto elemSol = elementSolution(element, sol, gg);
        return evalSolution(element, element.geometry(), gg, elemSol, scvf.ipGlobal())[pressureIdx];
    }

    /*!
     * \brief Returns the pressure at a given sub control volume
     */
    Scalar pressure(const Element<freeFlowMomentumIndex>& element,
                    const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                    const SubControlVolume<freeFlowMomentumIndex>& scv,
                    const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !this->isTransient_));
        const auto& gg = this->problem(freeFlowMassIndex).gridGeometry();
        const auto& sol = considerPreviousTimeStep ? (*prevSol_)[freeFlowMassIndex]
                                                   :  this->curSol(freeFlowMassIndex);
        const auto elemSol = elementSolution(element, sol, gg);
        return evalSolution(element, element.geometry(), gg, elemSol, scv.dofPosition())[pressureIdx];
    }

    /*!
     * \brief Returns the density at a given sub control volume face.
     */
    Scalar density(const Element<freeFlowMomentumIndex>& element,
                   const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                   const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                   const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !this->isTransient_));
        bindCouplingContext_(Dune::index_constant<freeFlowMomentumIndex>(), element, fvGeometry.elementIndex());

        if constexpr (MassDiscretizationMethod{} == DiscretizationMethods::cctpfa)
        {
            const auto eIdx = fvGeometry.elementIndex();
            const auto& scv = this->momentumCouplingContext_()[0].fvGeometry.scv(eIdx);

            const auto& volVars = considerPreviousTimeStep ?
                this->momentumCouplingContext_()[0].prevElemVolVars[scv]
                : this->momentumCouplingContext_()[0].curElemVolVars[scv];

            return volVars.density();
        }
        else if constexpr (MassDiscretizationMethod{} == DiscretizationMethods::box
                           || MassDiscretizationMethod{} == DiscretizationMethods::fcdiamond)
        {
            // TODO: cache the shape values when Box method is used
            using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
            const auto& localBasis = this->momentumCouplingContext_()[0].fvGeometry.feLocalBasis();
            std::vector<ShapeValue> shapeValues;
            const auto ipLocal = element.geometry().local(scvf.ipGlobal());
            localBasis.evaluateFunction(ipLocal, shapeValues);

            Scalar rho = 0.0;
            for (const auto& scv : scvs(this->momentumCouplingContext_()[0].fvGeometry))
            {
                const auto& volVars = considerPreviousTimeStep ?
                    this->momentumCouplingContext_()[0].prevElemVolVars[scv]
                    : this->momentumCouplingContext_()[0].curElemVolVars[scv];
                rho += volVars.density()*shapeValues[scv.indexInElement()][0];
            }

            return rho;
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                "Density interpolation for discretization scheme " << MassDiscretizationMethod{}
            );
    }

    /*!
     * \brief Returns the density at a given sub control volume.
     */
    Scalar density(const Element<freeFlowMomentumIndex>& element,
                   const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                   const SubControlVolume<freeFlowMomentumIndex>& scv,
                   const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !this->isTransient_));
        bindCouplingContext_(Dune::index_constant<freeFlowMomentumIndex>(), element, scv.elementIndex());

        if constexpr (MassDiscretizationMethod{} == DiscretizationMethods::cctpfa)
        {
            const auto eIdx = scv.elementIndex();
            const auto& scvI = this->momentumCouplingContext_()[0].fvGeometry.scv(eIdx);

            const auto& volVars = considerPreviousTimeStep ?
                this->momentumCouplingContext_()[0].prevElemVolVars[scvI]
                : this->momentumCouplingContext_()[0].curElemVolVars[scvI];

            return volVars.density();
        }
        else if constexpr (MassDiscretizationMethod{} == DiscretizationMethods::box
                           || MassDiscretizationMethod{} == DiscretizationMethods::fcdiamond)
        {
            // TODO: cache the shape values when Box method is used
            using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
            const auto& localBasis = this->momentumCouplingContext_()[0].fvGeometry.feLocalBasis();
            std::vector<ShapeValue> shapeValues;
            const auto ipLocal = element.geometry().local(scv.dofPosition());
            localBasis.evaluateFunction(ipLocal, shapeValues);

            Scalar rho = 0.0;
            for (const auto& scvI : scvs(this->momentumCouplingContext_()[0].fvGeometry))
            {
                const auto& volVars = considerPreviousTimeStep ?
                    this->momentumCouplingContext_()[0].prevElemVolVars[scvI]
                    : this->momentumCouplingContext_()[0].curElemVolVars[scvI];
                rho += volVars.density()*shapeValues[scvI.indexInElement()][0];
            }
            return rho;
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                "Density interpolation for discretization scheme " << MassDiscretizationMethod{}
            );
    }

    /*!
     * \brief Returns the effective viscosity at a given sub control volume face.
     */
    Scalar effectiveViscosity(const Element<freeFlowMomentumIndex>& element,
                              const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                              const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                              const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !this->isTransient_));
        bindCouplingContext_(Dune::index_constant<freeFlowMomentumIndex>(), element, fvGeometry.elementIndex());

        if constexpr (MassDiscretizationMethod{} == DiscretizationMethods::cctpfa)
        {
            const auto eIdx = fvGeometry.elementIndex();
            const auto& scv = this->momentumCouplingContext_()[0].fvGeometry.scv(eIdx);
            const auto& volVars = considerPreviousTimeStep ?
                this->momentumCouplingContext_()[0].prevElemVolVars[scv]
                : this->momentumCouplingContext_()[0].curElemVolVars[scv];
            return volVars.viscosity();
        }
        else if constexpr (MassDiscretizationMethod{} == DiscretizationMethods::box
                           || MassDiscretizationMethod{} == DiscretizationMethods::fcdiamond)
        {
            // TODO: cache the shape values when Box method is used
            using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
            const auto& localBasis = this->momentumCouplingContext_()[0].fvGeometry.feLocalBasis();
            std::vector<ShapeValue> shapeValues;
            const auto ipLocal = element.geometry().local(scvf.ipGlobal());
            localBasis.evaluateFunction(ipLocal, shapeValues);

            Scalar mu = 0.0;
            for (const auto& scv : scvs(this->momentumCouplingContext_()[0].fvGeometry))
            {
                const auto& volVars = considerPreviousTimeStep ?
                    this->momentumCouplingContext_()[0].prevElemVolVars[scv]
                    : this->momentumCouplingContext_()[0].curElemVolVars[scv];
                mu += volVars.viscosity()*shapeValues[scv.indexInElement()][0];
            }

            return mu;
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                "Viscosity interpolation for discretization scheme " << MassDiscretizationMethod{}
            );
    }

    /*!
     * \brief Returns the effective viscosity at a given sub control volume.
     */
    Scalar effectiveViscosity(const Element<freeFlowMomentumIndex>& element,
                              const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                              const SubControlVolume<freeFlowMomentumIndex>& scv,
                              const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !this->isTransient_));
        bindCouplingContext_(Dune::index_constant<freeFlowMomentumIndex>(), element, fvGeometry.elementIndex());

        if constexpr (MassDiscretizationMethod{} == DiscretizationMethods::cctpfa)
        {
            const auto eIdx = fvGeometry.elementIndex();
            const auto& scvI = this->momentumCouplingContext_()[0].fvGeometry.scv(eIdx);
            const auto& volVars = considerPreviousTimeStep ?
                this->momentumCouplingContext_()[0].prevElemVolVars[scvI]
                : this->momentumCouplingContext_()[0].curElemVolVars[scvI];
            return volVars.viscosity();
        }
        else if constexpr (MassDiscretizationMethod{} == DiscretizationMethods::box
                           || MassDiscretizationMethod{} == DiscretizationMethods::fcdiamond)
        {
            // TODO: cache the shape values when Box method is used
            using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
            const auto& localBasis = this->momentumCouplingContext_()[0].fvGeometry.feLocalBasis();
            std::vector<ShapeValue> shapeValues;
            const auto ipLocal = element.geometry().local(scv.dofPosition());
            localBasis.evaluateFunction(ipLocal, shapeValues);

            Scalar mu = 0.0;
            for (const auto& scvI : scvs(this->momentumCouplingContext_()[0].fvGeometry))
            {
                const auto& volVars = considerPreviousTimeStep ?
                    this->momentumCouplingContext_()[0].prevElemVolVars[scvI]
                    : this->momentumCouplingContext_()[0].curElemVolVars[scvI];
                mu += volVars.viscosity()*shapeValues[scvI.indexInElement()][0];
            }

            return mu;
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                "Viscosity interpolation for discretization scheme " << MassDiscretizationMethod{}
            );
    }

     /*!
     * \brief Returns the velocity at a given sub control volume face.
     */
    VelocityVector faceVelocity(const Element<freeFlowMassIndex>& element,
                                const SubControlVolumeFace<freeFlowMassIndex>& scvf) const
    {
        // TODO: optimize this function for tpfa where the scvf ip coincides with the dof location

        const auto eIdx = this->problem(freeFlowMassIndex).gridGeometry().elementMapper().index(element);
        bindCouplingContext_(Dune::index_constant<freeFlowMassIndex>(), element, eIdx);

        const auto& fvGeometry = this->massAndEnergyCouplingContext_()[0].fvGeometry;
        const auto& localBasis = fvGeometry.feLocalBasis();

        std::vector<ShapeValue> shapeValues;
        const auto ipLocal = element.geometry().local(scvf.ipGlobal());
        localBasis.evaluateFunction(ipLocal, shapeValues);

        // interpolate velocity at scvf
        VelocityVector velocity(0.0);
        for (const auto& scv : scvs(fvGeometry))
            velocity.axpy(shapeValues[scv.localDofIndex()][0], this->curSol(freeFlowMomentumIndex)[scv.dofIndex()]);

        return velocity;
    }

    /*!
     * \brief Returns the velocity at the element center.
     */
    VelocityVector elementVelocity(const FVElementGeometry<freeFlowMassIndex>& fvGeometry) const
    {
        bindCouplingContext_(Dune::index_constant<freeFlowMassIndex>(), fvGeometry.element());

        const auto& momentumFvGeometry = this->massAndEnergyCouplingContext_()[0].fvGeometry;
        const auto& localBasis = momentumFvGeometry.feLocalBasis();

        // interpolate velocity at scvf
        VelocityVector velocity(0.0);
        std::vector<ShapeValue> shapeValues;
        localBasis.evaluateFunction(referenceElement(fvGeometry.element()).position(0,0), shapeValues);

        for (const auto& scv : scvs(momentumFvGeometry))
            velocity.axpy(shapeValues[scv.localDofIndex()][0], this->curSol(freeFlowMomentumIndex)[scv.dofIndex()]);

        return velocity;
    }

    /*!
     * \brief The coupling stencil of domain I, i.e. which domain J DOFs
     *        the given domain I element's residual depends on.
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
     *        that couple with / influence the element residual of the given element of domain i
     *
     * \param domainI the domain index of domain i
     * \param elementI the coupled element of domain í
     * \param domainJ the domain index of domain j
     */
    const CouplingStencilType& couplingStencil(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                               const Element<freeFlowMomentumIndex>& elementI,
                                               Dune::index_constant<freeFlowMassIndex> domainJ) const
    {
        const auto eIdx = this->problem(freeFlowMomentumIndex).gridGeometry().elementMapper().index(elementI);
        return momentumToMassAndEnergyStencils_[eIdx];
    }

    // \}

    /*!
     * \name member functions concerning variable caching for element residual evaluations
     */
    // \{

    /*!
     * \ingroup MultiDomain
     * \brief updates all data and variables that are necessary to evaluate the residual of the element of domain i
     *        this is called whenever one of the primary variables that the element residual depends on changes in domain j
     *
     * \param domainI the domain index of domain i
     * \param localAssemblerI the local assembler assembling the element residual of an element of domain i
     * \param domainJ the domain index of domain j
     * \param dofIdxGlobalJ the index of the degree of freedom of domain j whose solution changed
     * \param priVarsJ the new solution at the degree of freedom of domain j with index dofIdxGlobalJ
     * \param pvIdxJ the index of the primary variable of domain j which has been updated
     *
     * \note this concerns all data that is used in the evaluation of the element residual and depends on
     *       the primary variables at the degree of freedom location with index dofIdxGlobalJ
     * \note  the element whose residual is to be evaluated can be retrieved from the local assembler
     *        as localAssemblerI.element()
     * \note  per default, we update the solution vector, if the element residual of domain i depends on more than
     *        the primary variables of domain j update the other dependent data here by overloading this function
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

        if constexpr (MassDiscretizationMethod{} == DiscretizationMethods::cctpfa)
        {
            if constexpr (domainI == freeFlowMomentumIndex && domainJ == freeFlowMassIndex)
            {
                bindCouplingContext_(domainI, localAssemblerI.element());

                const auto& problem = this->problem(domainJ);
                const auto& deflectedElement = problem.gridGeometry().element(dofIdxGlobalJ);
                const auto elemSol = elementSolution(deflectedElement, this->curSol(domainJ), problem.gridGeometry());
                const auto& fvGeometry = momentumCouplingContext_()[0].fvGeometry;
                const auto& scv = fvGeometry.scv(dofIdxGlobalJ);

                if constexpr (ElementVolumeVariables<freeFlowMassIndex>::GridVolumeVariables::cachingEnabled)
                    gridVars_(freeFlowMassIndex).curGridVolVars().volVars(scv).update(std::move(elemSol), problem, deflectedElement, scv);
                else
                    momentumCouplingContext_()[0].curElemVolVars[scv].update(std::move(elemSol), problem, deflectedElement, scv);
            }
        }
        else if constexpr (MassDiscretizationMethod{} == DiscretizationMethods::box
                           || MassDiscretizationMethod{} == DiscretizationMethods::fcdiamond)
        {
            if constexpr (domainI == freeFlowMomentumIndex && domainJ == freeFlowMassIndex)
            {
                bindCouplingContext_(domainI, localAssemblerI.element());

                const auto& problem = this->problem(domainJ);
                const auto& deflectedElement = problem.gridGeometry().element(this->momentumCouplingContext_()[0].eIdx);
                const auto elemSol = elementSolution(deflectedElement, this->curSol(domainJ), problem.gridGeometry());
                const auto& fvGeometry = this->momentumCouplingContext_()[0].fvGeometry;

                for (const auto& scv : scvs(fvGeometry))
                {
                    if(scv.dofIndex() == dofIdxGlobalJ)
                    {
                        if constexpr (ElementVolumeVariables<freeFlowMassIndex>::GridVolumeVariables::cachingEnabled)
                            this->gridVars_(freeFlowMassIndex).curGridVolVars().volVars(scv).update(std::move(elemSol), problem, deflectedElement, scv);
                        else
                           this->momentumCouplingContext_()[0].curElemVolVars[scv].update(std::move(elemSol), problem, deflectedElement, scv);
                    }
                }
            }
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                "Context update for discretization scheme " << MassDiscretizationMethod{}
            );
    }

    // \}

    /*!
     * \brief Compute colors for multithreaded assembly
     */
    void computeColorsForAssembly()
    {
        // use coloring of the momentum discretization for both domains
        elementSets_ = computeColoring(this->problem(freeFlowMomentumIndex).gridGeometry()).sets;
    }

    /*!
     * \brief Execute assembly kernel in parallel
     *
     * \param domainId the domain index of domain i
     * \param assembleElement kernel function to execute for one element
     */
    template<std::size_t i, class AssembleElementFunc>
    void assembleMultithreaded(Dune::index_constant<i> domainId, AssembleElementFunc&& assembleElement) const
    {
        if (elementSets_.empty())
            DUNE_THROW(Dune::InvalidStateException, "Call computeColorsForAssembly before assembling in parallel!");

        // make this element loop run in parallel
        // for this we have to color the elements so that we don't get
        // race conditions when writing into the global matrix
        // each color can be assembled using multiple threads
        const auto& grid = this->problem(freeFlowMassIndex).gridGeometry().gridView().grid();
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
    void bindCouplingContext_(Dune::index_constant<freeFlowMomentumIndex> domainI,
                              const Element<freeFlowMomentumIndex>& elementI) const
    {
        // The call to this->problem() is expensive because of std::weak_ptr (see base class). Here we try to avoid it if possible.
        if (momentumCouplingContext_().empty())
            bindCouplingContext_(domainI, elementI, this->problem(freeFlowMomentumIndex).gridGeometry().elementMapper().index(elementI));
        else
            bindCouplingContext_(domainI, elementI, momentumCouplingContext_()[0].fvGeometry.gridGeometry().elementMapper().index(elementI));
    }

    void bindCouplingContext_(Dune::index_constant<freeFlowMomentumIndex> domainI,
                              const Element<freeFlowMomentumIndex>& elementI,
                              const std::size_t eIdx) const
    {
        if (momentumCouplingContext_().empty())
        {
            auto fvGeometry = localView(this->problem(freeFlowMassIndex).gridGeometry());
            fvGeometry.bind(elementI);

            auto curElemVolVars = localView(gridVars_(freeFlowMassIndex).curGridVolVars());
            curElemVolVars.bind(elementI, fvGeometry, this->curSol(freeFlowMassIndex));

            auto prevElemVolVars = isTransient_ ? localView(gridVars_(freeFlowMassIndex).prevGridVolVars())
                                                : localView(gridVars_(freeFlowMassIndex).curGridVolVars());

            if (isTransient_)
                prevElemVolVars.bindElement(elementI, fvGeometry, (*prevSol_)[freeFlowMassIndex]);

            momentumCouplingContext_().emplace_back(MomentumCouplingContext{std::move(fvGeometry), std::move(curElemVolVars), std::move(prevElemVolVars), eIdx});
        }
        else if (eIdx != momentumCouplingContext_()[0].eIdx)
        {
            momentumCouplingContext_()[0].eIdx = eIdx;
            momentumCouplingContext_()[0].fvGeometry.bind(elementI);
            momentumCouplingContext_()[0].curElemVolVars.bind(elementI, momentumCouplingContext_()[0].fvGeometry, this->curSol(freeFlowMassIndex));

            if (isTransient_)
                momentumCouplingContext_()[0].prevElemVolVars.bindElement(elementI, momentumCouplingContext_()[0].fvGeometry, (*prevSol_)[freeFlowMassIndex]);
        }
    }

    void bindCouplingContext_(Dune::index_constant<freeFlowMassIndex> domainI,
                              const Element<freeFlowMassIndex>& elementI) const
    {
        // The call to this->problem() is expensive because of std::weak_ptr (see base class). Here we try to avoid it if possible.
        if (massAndEnergyCouplingContext_().empty())
            bindCouplingContext_(domainI, elementI, this->problem(freeFlowMassIndex).gridGeometry().elementMapper().index(elementI));
        else
            bindCouplingContext_(domainI, elementI, massAndEnergyCouplingContext_()[0].fvGeometry.gridGeometry().elementMapper().index(elementI));
    }

    void bindCouplingContext_(Dune::index_constant<freeFlowMassIndex> domainI,
                              const Element<freeFlowMassIndex>& elementI,
                              const std::size_t eIdx) const
    {
        if (massAndEnergyCouplingContext_().empty())
        {
            const auto& gridGeometry = this->problem(freeFlowMomentumIndex).gridGeometry();
            auto fvGeometry = localView(gridGeometry);
            fvGeometry.bindElement(elementI);
            massAndEnergyCouplingContext_().emplace_back(std::move(fvGeometry), eIdx);
        }
        else if (eIdx != massAndEnergyCouplingContext_()[0].eIdx)
        {
            massAndEnergyCouplingContext_()[0].eIdx = eIdx;
            massAndEnergyCouplingContext_()[0].fvGeometry.bindElement(elementI);
        }
    }

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
        const auto& momentumGridGeometry = this->problem(freeFlowMomentumIndex).gridGeometry();
        const auto& massGridGeometry = this->problem(freeFlowMassIndex).gridGeometry();
        auto momentumFvGeometry = localView(momentumGridGeometry);
        auto massFvGeometry = localView(massGridGeometry);

        massAndEnergyToMomentumStencils_.clear();
        massAndEnergyToMomentumStencils_.resize(massGridGeometry.gridView().size(0));

        momentumToMassAndEnergyStencils_.clear();
        momentumToMassAndEnergyStencils_.resize(momentumGridGeometry.gridView().size(0));

        assert(massAndEnergyToMomentumStencils_.size() == momentumToMassAndEnergyStencils_.size());

        for (const auto& element : elements(momentumGridGeometry.gridView()))
        {
            momentumFvGeometry.bindElement(element);
            massFvGeometry.bindElement(element);
            const auto eIdx = momentumFvGeometry.elementIndex();

            for (const auto& scv : scvs(momentumFvGeometry))
                massAndEnergyToMomentumStencils_[eIdx].push_back(scv.dofIndex());

            for (const auto& scv : scvs(massFvGeometry))
                momentumToMassAndEnergyStencils_[eIdx].push_back(scv.dofIndex());
        }
    }

    CouplingStencilType emptyStencil_;
    std::vector<CouplingStencilType> momentumToMassAndEnergyStencils_;
    std::vector<CouplingStencilType> massAndEnergyToMomentumStencils_;

    // the coupling context exists for each thread
    // TODO this is a bad pattern, just like mutable caches
    // we should really construct and pass the context and not store it globally
    std::vector<MomentumCouplingContext>& momentumCouplingContext_() const
    {
        thread_local static std::vector<MomentumCouplingContext> c;
        return c;
    }

    // the coupling context exists for each thread
    std::vector<MassAndEnergyCouplingContext>& massAndEnergyCouplingContext_() const
    {
        thread_local static std::vector<MassAndEnergyCouplingContext> c;
        return c;
    }

    //! A tuple of std::shared_ptrs to the grid variables of the sub problems
    GridVariablesTuple gridVariables_;

    const SolutionVector* prevSol_;
    bool isTransient_;

    std::deque<std::vector<ElementSeed<freeFlowMomentumIndex>>> elementSets_;
};

//! TODO The infrastructure for multithreaded assembly is implemented (see code in the class above) but the current implementation seems to have a bug and may cause race conditions. The result is different when running in parallel. After this has been fixed activate multithreaded assembly by inheriting from std::true_type here.
template<class T>
struct CouplingManagerSupportsMultithreadedAssembly<PQ1BubbleFreeFlowCouplingManager<T>>
: public std::false_type {};

} // end namespace Dumux

#endif
