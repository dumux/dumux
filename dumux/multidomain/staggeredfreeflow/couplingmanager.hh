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
 * \ingroup MultiDomain
 * \brief The interface of the coupling manager for multi domain problems
 */

#ifndef DUMUX_MULTIDOMAIN_STAGGERED_FREEFLOW_COUPLING_MANAGER_HH
#define DUMUX_MULTIDOMAIN_STAGGERED_FREEFLOW_COUPLING_MANAGER_HH

#include <memory>
#include <tuple>
#include <vector>
#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/typetraits.hh>

#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup MultiDomain
 * \brief The interface of the coupling manager for multi domain problems
 */
template<class Traits>
class StaggeredFreeFlowCouplingManagerBase : public CouplingManager<Traits>
{
public:
    static constexpr auto freeFlowMomentumIdx = typename Traits::template SubDomain<0>::Index();
    static constexpr auto freeFlowMassIdx = typename Traits::template SubDomain<1>::Index();
protected:
    template<std::size_t id> using SubDomainTypeTag = typename Traits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
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

    using ParentType = CouplingManager<Traits>;

    using CouplingStencilType = std::vector<std::size_t>;

    using GridVariablesTuple = typename Traits::template TupleOfSharedPtr<GridVariables>;

    using FluidSystem = typename VolumeVariables<freeFlowMassIdx>::FluidSystem;
    // static_assert(std::is_same_v<FluidSystem, typename VolumeVariables<freeFlowMomentumIdx>::FluidSystem>);

    using VelocityVector = typename SubControlVolumeFace<freeFlowMassIdx>::GlobalPosition;
    static_assert(std::is_same_v<VelocityVector, typename SubControlVolumeFace<freeFlowMomentumIdx>::GlobalPosition>);
    struct MomentumCouplingContext
    {
        FVElementGeometry<freeFlowMassIdx> fvGeometry;
        ElementVolumeVariables<freeFlowMassIdx> curElemVolVars;
        ElementVolumeVariables<freeFlowMassIdx> prevElemVolVars;
        std::size_t eIdx;
    };

    struct MassAndEnergyCouplingContext
    {
        FVElementGeometry<freeFlowMomentumIdx> fvGeometry;
        std::size_t eIdx;
    };

public:

    static constexpr auto pressureIdx = VolumeVariables<freeFlowMassIdx>::Indices::pressureIdx;

    /*!
     * \ingroup MultiDomain
     * \brief evaluates the element residual of a coupled element of domain i which depends on the variables
     *        at the degree of freedom with index dofIdxGlobalJ of domain j
     *
     * \param domainI the domain index of domain i
     * \param localAssemblerI the local assembler assembling the element residual of an element of domain i
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
    decltype(auto) evalCouplingResidual(Dune::index_constant<freeFlowMomentumIdx> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        const SubControlVolume<freeFlowMomentumIdx>& scvI,
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
            assert(isTransient_);
            localResidual.evalStorage(residual, problem, element, fvGeometry, prevElemVolVars, curElemVolVars, scvI);
        }

        return residual;
    }


    //! TODO: this is just a prototype. May be removed after some testing
    Scalar extrapolatedPressure(const Element<freeFlowMomentumIdx>& element,
                                const FVElementGeometry<freeFlowMomentumIdx>& fvGeometry,
                                const SubControlVolumeFace<freeFlowMomentumIdx>& scvf) const
    {
        const auto ownCellPressure = pressure(element, fvGeometry, scvf);
        assert (scvf.boundary() && scvf.isFrontal());

        bindCouplingContext(Dune::index_constant<freeFlowMomentumIdx>(), element, fvGeometry.elementIndex());

        for (const auto is : intersections(fvGeometry.gridGeometry().gridView(), element))
        {
            if ((is.centerUnitOuterNormal() * scvf.unitOuterNormal() + 1) < 1e-9)
            {
                const auto& outsideElement = is.outside();
                const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(outsideElement);
                const auto& outsideScv = momentumCouplingContext_[0].fvGeometry.scv(eIdx);
                const auto p = momentumCouplingContext_[0].curElemVolVars[outsideScv].pressure() - this->problem(freeFlowMassIdx).initial(outsideElement)[pressureIdx];
                const auto slope = (ownCellPressure - p) / (element.geometry().center() - outsideScv.center()).two_norm();

                return ownCellPressure + slope * (element.geometry().center() - scvf.center()).two_norm(); // only works if boundary is on the right
            }
        }
        DUNE_THROW(Dune::InvalidStateException, "No intersection found");
    }

    /*!
     * \name member functions concerning the coupling stencils
     */
    // \{

    /*!
     * \brief Returns the pressure at a given sub control volume face.
     */
    Scalar pressure(const Element<freeFlowMomentumIdx>& element,
                    const FVElementGeometry<freeFlowMomentumIdx>& fvGeometry,
                    const SubControlVolumeFace<freeFlowMomentumIdx>& scvf) const
    {
        return this->curSol()[freeFlowMassIdx][fvGeometry.elementIndex()][pressureIdx];
    }

    /*!
     * \brief Returns the pressure at the center of a sub control volume corresponding to a given sub control volume face.
     *        This is used for setting a Dirichlet pressure for the mass model when a fixed pressure for the momentum balance is set at another
     *        boundary. Since the the pressure at the given scvf is solution-dependent and thus unknown a priori, we just use the value
     *        of the interior cell here.
     */
    Scalar cellPressure(const Element<freeFlowMassIdx>& element,
                        const SubControlVolumeFace<freeFlowMassIdx>& scvf) const
    {
        return this->curSol()[freeFlowMassIdx][scvf.insideScvIdx()][pressureIdx];
    }

    /*!
     * \brief Returns the density at a given sub control volume face.
     */
    Scalar density(const Element<freeFlowMomentumIdx>& element,
                   const FVElementGeometry<freeFlowMomentumIdx>& fvGeometry,
                   const SubControlVolumeFace<freeFlowMomentumIdx>& scvf,
                   const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !isTransient_));
        bindCouplingContext(Dune::index_constant<freeFlowMomentumIdx>(), element, fvGeometry.elementIndex());
        const auto& insideMomentumScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideMassScv = momentumCouplingContext_[0].fvGeometry.scv(insideMomentumScv.elementIndex());

        auto rho = [&](const auto& elemVolVars)
        {
            if (scvf.boundary())
                return elemVolVars[insideMassScv].density();
            else
            {
                const auto& outsideMomentumScv = fvGeometry.scv(scvf.outsideScvIdx());
                const auto& outsideMassScv = momentumCouplingContext_[0].fvGeometry.scv(outsideMomentumScv.elementIndex());
                // TODO distance weighting
                return 0.5*(elemVolVars[insideMassScv].density() + elemVolVars[outsideMassScv].density());
            }
        };

        return considerPreviousTimeStep ? rho(momentumCouplingContext_[0].prevElemVolVars)
                                        : rho(momentumCouplingContext_[0].curElemVolVars);
    }

    auto getInsideAndOutsideDensity(const Element<freeFlowMomentumIdx>& element,
                                    const FVElementGeometry<freeFlowMomentumIdx>& fvGeometry,
                                    const SubControlVolumeFace<freeFlowMomentumIdx>& scvf,
                                    const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !isTransient_));
        bindCouplingContext(Dune::index_constant<freeFlowMomentumIdx>(), element, fvGeometry.elementIndex());
        const auto& insideMomentumScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideMassScv = momentumCouplingContext_[0].fvGeometry.scv(insideMomentumScv.elementIndex());

        auto result = [&](const auto& elemVolVars)
        {
            if (scvf.boundary())
                return std::make_pair(elemVolVars[insideMassScv].density(), elemVolVars[insideMassScv].density());
            else
            {
                const auto& outsideMomentumScv = fvGeometry.scv(scvf.outsideScvIdx());
                const auto& outsideMassScv = momentumCouplingContext_[0].fvGeometry.scv(outsideMomentumScv.elementIndex());
                return std::make_pair(elemVolVars[insideMassScv].density(), elemVolVars[outsideMassScv].density());
            }
        };

        return considerPreviousTimeStep ? result(momentumCouplingContext_[0].prevElemVolVars)
                                        : result(momentumCouplingContext_[0].curElemVolVars);
    }

    /*!
     * \brief Returns the density at a given sub control volume.
     */
    Scalar density(const Element<freeFlowMomentumIdx>& element,
                   const SubControlVolume<freeFlowMomentumIdx>& scv,
                   const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !isTransient_));
        bindCouplingContext(Dune::index_constant<freeFlowMomentumIdx>(), element, scv.elementIndex());
        const auto& massScv = (*scvs(momentumCouplingContext_[0].fvGeometry).begin());

        return considerPreviousTimeStep ? momentumCouplingContext_[0].prevElemVolVars[massScv].density()
                                        : momentumCouplingContext_[0].curElemVolVars[massScv].density();
    }

    /*!
     * \brief Returns the pressure at a given sub control volume face.
     */
    Scalar effectiveViscosity(const Element<freeFlowMomentumIdx>& element,
                              const FVElementGeometry<freeFlowMomentumIdx>& fvGeometry,
                              const SubControlVolumeFace<freeFlowMomentumIdx>& scvf) const
    {
        bindCouplingContext(Dune::index_constant<freeFlowMomentumIdx>(), element, fvGeometry.elementIndex());

        const auto& insideMomentumScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideMassScv = momentumCouplingContext_[0].fvGeometry.scv(insideMomentumScv.elementIndex());

        if (scvf.boundary())
            return momentumCouplingContext_[0].curElemVolVars[insideMassScv].viscosity();

        const auto& outsideMomentumScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& outsideMassScv = momentumCouplingContext_[0].fvGeometry.scv(outsideMomentumScv.elementIndex());

        auto mu = [&](const auto& elemVolVars)
        {
            // TODO distance weighting
            return 0.5*(elemVolVars[insideMassScv].viscosity() + elemVolVars[outsideMassScv].viscosity());
        };

        return mu(momentumCouplingContext_[0].curElemVolVars);
    }

     /*!
     * \brief Returns the velocity at a given sub control volume face.
     */
    virtual VelocityVector faceVelocity(const Element<freeFlowMassIdx>& element,
                                const SubControlVolumeFace<freeFlowMassIdx>& scvf) const = 0;

    /*!
     * \brief The coupling stencil of domain I, i.e. which domain J DOFs
     *        the given domain I element's residual depends on.
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencilType& couplingStencil(Dune::index_constant<i> domainI,
                                               const Element<i>& element,
                                               Dune::index_constant<j> domainJ) const
    { return emptyStencil_; }

    /*!
     * \brief returns an iteratable container of all indices of degrees of freedom of domain j
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
    const CouplingStencilType couplingStencil(Dune::index_constant<freeFlowMassIdx> domainI,
                                              const Element<freeFlowMassIdx>& elementI,
                                              Dune::index_constant<freeFlowMomentumIdx> domainJ) const
    {
        const auto eIdx = this->problem(freeFlowMassIdx).gridGeometry().elementMapper().index(elementI);
        return massAndEnergyToMomentumStencils_[eIdx];
    }

    /*!
     * \brief returns an iteratable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the residual of the given sub-control volume of domain i
     *
     * \param domainI the domain index of domain i
     * \param elementI the coupled element of domain í
     * \param scvI the sub-control volume of domain i
     * \param domainJ the domain index of domain j
     */
    CouplingStencilType couplingStencil(Dune::index_constant<freeFlowMomentumIdx> domainI,
                                        const Element<freeFlowMomentumIdx>& elementI,
                                        const SubControlVolume<freeFlowMomentumIdx>& scvI,
                                        Dune::index_constant<freeFlowMassIdx> domainJ) const
    {
        return momentumToMassAndEnergyStencils_[scvI.index()];
    }

    /*!
     * \brief extend the jacobian pattern of the diagonal block of domain i
     *        by those entries that are not already in the uncoupled pattern
     * \note per default we do not add such additional dependencies
     * \note Such additional dependencies can arise from the coupling, e.g. if a coupling source
     *       term depends on a non-local average of a quantity of the same domain
     * \warning if you overload this also implement evalAdditionalDomainDerivatives
     */
    template<std::size_t id, class JacobianPattern>
    void extendJacobianPattern(Dune::index_constant<id> domainI, JacobianPattern& pattern) const
    {}

    // \}

    /*!
     * \name member functions concerning variable caching for element residual evaluations
     */
    // \{

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of the element of domain i
     *
     * \param domainI the domain index of domain i
     * \param elementI the element whose residual we are assemling next
     * \param assembler the multidomain assembler for access to all data necessary for the assembly of all domains
     *
     * \note this concerns all data that is used in the evaluation of the element residual and depends on one of the
     *       degrees of freedom returned by CouplingManager::couplingStencil
     * \note every coupled element residual depends at least on the solution of another domain, that why we always store a
     *       copy of the solution vector in the coupling manager, hence, in case the element residual
     *       only depends on primary variables of the other domain this function does nothing
     * \note overload this function in case the element residual depends on more than the primary variables of domain j
     */
    template<std::size_t i, class Assembler>
    void bindCouplingContext(Dune::index_constant<i> domainI,
                             const Element<i>& elementI,
                             const Assembler& assembler)
    {}

    void bindCouplingContext(Dune::index_constant<freeFlowMomentumIdx> domainI,
                             const Element<freeFlowMomentumIdx>& elementI) const
    {
        // The call to this->problem() is expensive because of std::weak_ptr (see base class). Here we try to avoid it if possible.
        if (momentumCouplingContext_.empty())
            bindCouplingContext(domainI, elementI, this->problem(freeFlowMomentumIdx).gridGeometry().elementMapper().index(elementI));
        else
            bindCouplingContext(domainI, elementI, momentumCouplingContext_[0].fvGeometry.gridGeometry().elementMapper().index(elementI));
    }

    void bindCouplingContext(Dune::index_constant<freeFlowMomentumIdx> domainI,
                             const Element<freeFlowMomentumIdx>& elementI,
                             const std::size_t eIdx) const
    {
        if (momentumCouplingContext_.empty())
        {
            auto fvGeometry = localView(this->problem(freeFlowMassIdx).gridGeometry());
            fvGeometry.bind(elementI);

            auto curElemVolVars = localView(gridVars_(freeFlowMassIdx).curGridVolVars());
            curElemVolVars.bind(elementI, fvGeometry, this->curSol()[freeFlowMassIdx]);

            auto prevElemVolVars = isTransient_ ? localView(gridVars_(freeFlowMassIdx).prevGridVolVars())
                                                : localView(gridVars_(freeFlowMassIdx).curGridVolVars());

            if (isTransient_)
                prevElemVolVars.bindElement(elementI, fvGeometry, (*prevSol_)[freeFlowMassIdx]);

            momentumCouplingContext_.emplace_back(MomentumCouplingContext{std::move(fvGeometry), std::move(curElemVolVars), std::move(prevElemVolVars), eIdx});
        }
        else if (eIdx != momentumCouplingContext_[0].eIdx)
        {
            momentumCouplingContext_[0].eIdx = eIdx;
            momentumCouplingContext_[0].fvGeometry.bind(elementI);
            momentumCouplingContext_[0].curElemVolVars.bind(elementI, momentumCouplingContext_[0].fvGeometry, this->curSol()[freeFlowMassIdx]);

            if (isTransient_)
                momentumCouplingContext_[0].prevElemVolVars.bindElement(elementI, momentumCouplingContext_[0].fvGeometry, (*prevSol_)[freeFlowMassIdx]);
        }
    }

    void bindCouplingContext(Dune::index_constant<freeFlowMassIdx> domainI,
                             const Element<freeFlowMassIdx>& elementI) const
    {
        // The call to this->problem() is expensive because of std::weak_ptr (see base class). Here we try to avoid it if possible.
        if (massAndEnergyCouplingContext_.empty())
            bindCouplingContext(domainI, elementI, this->problem(freeFlowMassIdx).gridGeometry().elementMapper().index(elementI));
        else
            bindCouplingContext(domainI, elementI, massAndEnergyCouplingContext_[0].fvGeometry.gridGeometry().elementMapper().index(elementI));
    }

    void bindCouplingContext(Dune::index_constant<freeFlowMassIdx> domainI,
                             const Element<freeFlowMassIdx>& elementI,
                             const std::size_t eIdx) const
    {
        if (massAndEnergyCouplingContext_.empty())
        {
            auto fvGeometry = localView(this->problem(freeFlowMomentumIdx).gridGeometry());
            fvGeometry.bindElement(elementI);

            massAndEnergyCouplingContext_.emplace_back(MassAndEnergyCouplingContext{std::move(fvGeometry), eIdx});
        }
        else if (eIdx != massAndEnergyCouplingContext_[0].eIdx)
        {
            massAndEnergyCouplingContext_[0].eIdx = eIdx;
            massAndEnergyCouplingContext_[0].fvGeometry.bindElement(elementI);
        }
    }

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
     * \note  per default, we udpate the solution vector, if the element residual of domain i depends on more than
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
        this->curSol()[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        if constexpr (domainI == freeFlowMomentumIdx && domainJ == freeFlowMassIdx)
        {
            bindCouplingContext(domainI, localAssemblerI.element());

            const auto& problem = this->problem(domainJ);
            const auto& deflectedElement = problem.gridGeometry().element(dofIdxGlobalJ);
            const auto elemSol = elementSolution(deflectedElement, this->curSol()[domainJ], problem.gridGeometry());
            const auto& fvGeometry = momentumCouplingContext_[0].fvGeometry;
            const auto scvIdxJ = dofIdxGlobalJ;
            const auto& scv = fvGeometry.scv(scvIdxJ);

            if constexpr (ElementVolumeVariables<freeFlowMassIdx>::GridVolumeVariables::cachingEnabled)
                gridVars_(freeFlowMassIdx).curGridVolVars().volVars(scv).update(std::move(elemSol), problem, deflectedElement, scv);
            else
                momentumCouplingContext_[0].curElemVolVars[scv].update(std::move(elemSol), problem, deflectedElement, scv);
        }
    }

    // \}

    /*!
     * \ingroup MultiDomain
     * \brief evaluates the element residual of a coupled element of domain i which depends on the variables
     *        at the degree of freedom with index dofIdxGlobalJ of domain j
     *
     * \param domainI the domain index of domain i
     * \param localAssemblerI the local assembler assembling the element residual of an element of domain i
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
    template<std::size_t i, std::size_t j, class LocalAssemblerI>
    decltype(auto) evalCouplingResidual(Dune::index_constant<i> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ) const
    {
        return localAssemblerI.evalLocalResidual();
    }

    /*!
     * \brief evaluate additional derivatives of the element residual of a domain with respect
     *        to dofs in the same domain that are not in the regular stencil (see CouplingManager::extendJacobianPattern)
     * \note Such additional dependencies can arise from the coupling, e.g. if a coupling source
     *       term depends on a non-local average of a quantity of the same domain
     */
    template<std::size_t i, class LocalAssemblerI, class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDomainDerivatives(Dune::index_constant<i> domainI,
                                         const LocalAssemblerI& localAssemblerI,
                                         const typename LocalAssemblerI::LocalResidual::ElementResidualVector& origResiduals,
                                         JacobianMatrixDiagBlock& A,
                                         GridVariables& gridVariables)
    {}

    /*!
     * \brief return the numeric epsilon used for deflecting primary variables of coupled domain i
     */
    template<std::size_t i>
    decltype(auto) numericEpsilon(Dune::index_constant<i>,
                                  const std::string& paramGroup) const
    {
        constexpr auto numEq = PrimaryVariables<i>::dimension;
        return NumericEpsilon<typename Traits::Scalar, numEq>(paramGroup);
    }

protected:

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

    CouplingStencilType emptyStencil_;
    std::vector<CouplingStencilType> momentumToMassAndEnergyStencils_;
    std::vector<CouplingStencilType> massAndEnergyToMomentumStencils_;
    mutable std::vector<MomentumCouplingContext> momentumCouplingContext_;
    mutable std::vector<MassAndEnergyCouplingContext> massAndEnergyCouplingContext_;

    /*!
    * \brief A tuple of std::shared_ptrs to the grid variables of the sub problems
    */
    GridVariablesTuple gridVariables_;

    const SolutionVector* prevSol_;
    bool isTransient_;
};

template<class Traits>
class StaggeredFreeFlowCouplingManager : public StaggeredFreeFlowCouplingManagerBase<Traits>
{
public:
    static constexpr auto freeFlowMomentumIdx = typename Traits::template SubDomain<0>::Index();
    static constexpr auto freeFlowMassIdx = typename Traits::template SubDomain<1>::Index();

private:
    template<std::size_t id> using SubDomainTypeTag = typename Traits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
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

    using FluidSystem = typename VolumeVariables<freeFlowMassIdx>::FluidSystem;

    using VelocityVector = typename SubControlVolumeFace<freeFlowMassIdx>::GlobalPosition;

    using ParentType = StaggeredFreeFlowCouplingManagerBase<Traits>;

public:
   /*!
     * \brief Methods to be accessed by main
     */
    // \{

    void init(std::shared_ptr<Problem<freeFlowMomentumIdx>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIdx>> massProblem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol)
    {
        this->setSubProblems(std::make_tuple(momentumProblem, massProblem));
        this->gridVariables_ = gridVariables;
        this->updateSolution(curSol);

        this->computeCouplingStencils_();
    }

    void init(std::shared_ptr<Problem<freeFlowMomentumIdx>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIdx>> massProblem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol,
              const SolutionVector& prevSol)
    {
        init(momentumProblem, massProblem, std::forward<GridVariablesTuple>(gridVariables), curSol);
        this->prevSol_ = &prevSol;
        this->isTransient_ = true;
    }

    // \}

    /*!
    * \brief Returns the velocity at a given sub control volume face.
    */
    VelocityVector faceVelocity(const Element<freeFlowMassIdx>& element,
                                const SubControlVolumeFace<freeFlowMassIdx>& scvf) const
    {
        this->bindCouplingContext(Dune::index_constant<freeFlowMassIdx>(), element, scvf.insideScvIdx()/*eIdx*/);
        const auto& scvJ = this->massAndEnergyCouplingContext_[0].fvGeometry.scv(scvf.index()/*corresponds to scvIdx of staggered*/);

        // create a unit normal vector oriented in positive coordinate direction
        auto velocity = scvf.unitOuterNormal();
        using std::abs;
        std::for_each(velocity.begin(), velocity.end(), [](auto& v){ v = abs(v); });

        // create the actual velocity vector
        velocity *= this->curSol()[freeFlowMomentumIdx][scvJ.dofIndex()];

        return velocity;
    }

private:
    void computeCouplingStencils_()
    {
        // TODO higher order
        const auto& momentumGridGeometry = this->problem(freeFlowMomentumIdx).gridGeometry();
        auto momentumFvGeometry = localView(momentumGridGeometry);
        this->massAndEnergyToMomentumStencils_.clear();
        this->massAndEnergyToMomentumStencils_.resize(momentumGridGeometry.gridView().size(0));

        this->momentumToMassAndEnergyStencils_.clear();
        this->momentumToMassAndEnergyStencils_.resize(momentumGridGeometry.numScv());

        for (const auto& element : elements(momentumGridGeometry.gridView()))
        {
            const auto eIdx = momentumGridGeometry.elementMapper().index(element);
            momentumFvGeometry.bindElement(element);
            for (const auto& scv : scvs(momentumFvGeometry))
            {
                this->massAndEnergyToMomentumStencils_[eIdx].push_back(scv.dofIndex());
                this->momentumToMassAndEnergyStencils_[scv.index()].push_back(eIdx);

                // extend the stencil for fluids with variable viscosity and density,
                if constexpr (FluidSystem::isCompressible(0/*phaseIdx*/))
                // if constexpr (FluidSystem::isCompressible(0/*phaseIdx*/) || !FluidSystem::viscosityIsConstant(0/*phaseIdx*/)) // TODO fix on master
                {
                    for (const auto& scvf : scvfs(momentumFvGeometry, scv))
                    {
                        if (scvf.isLateral() && !scvf.boundary())
                        {
                            const auto& outsideScv = momentumFvGeometry.scv(scvf.outsideScvIdx());
                            this->momentumToMassAndEnergyStencils_[scv.index()].push_back(outsideScv.elementIndex());
                        }
                    }
                }
            }
        }
    }

};

template<class Traits>
class DiamondFreeFlowCouplingManager : public StaggeredFreeFlowCouplingManagerBase<Traits>
{
public:
    static constexpr auto freeFlowMomentumIdx = typename Traits::template SubDomain<0>::Index();
    static constexpr auto freeFlowMassIdx = typename Traits::template SubDomain<1>::Index();

private:
    template<std::size_t id> using SubDomainTypeTag = typename Traits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
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

    using FluidSystem = typename VolumeVariables<freeFlowMassIdx>::FluidSystem;

    using VelocityVector = typename SubControlVolumeFace<freeFlowMassIdx>::GlobalPosition;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;

    using ParentType = StaggeredFreeFlowCouplingManagerBase<Traits>;

public:
   /*!
     * \brief Methods to be accessed by main
     */
    // \{

    void init(std::shared_ptr<Problem<freeFlowMomentumIdx>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIdx>> massProblem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol)
    {
        this->setSubProblems(std::make_tuple(momentumProblem, massProblem));
        this->gridVariables_ = gridVariables;
        this->updateSolution(curSol);

        this->computeCouplingStencils_();
    }

    void init(std::shared_ptr<Problem<freeFlowMomentumIdx>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIdx>> massProblem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol,
              const SolutionVector& prevSol)
    {
        init(momentumProblem, massProblem, std::forward<GridVariablesTuple>(gridVariables), curSol);
        this->prevSol_ = &prevSol;
        this->isTransient_ = true;
    }

    // \}

    /*!
    * \brief Returns the velocity at a given sub control volume face.
    */
    VelocityVector faceVelocity(const Element<freeFlowMassIdx>& element,
                                const SubControlVolumeFace<freeFlowMassIdx>& scvf) const
    {
        this->bindCouplingContext(Dune::index_constant<freeFlowMassIdx>(), element, scvf.insideScvIdx()/*eIdx*/);
        const auto& fvGeometry = this->massAndEnergyCouplingContext_[0];

        const auto& localBasis = fvGeometry.feLocalBasis();
        std::vector<ShapeValue> shapeValues;
        const auto ipLocal = element.geometry().local(scvf.ipGlobal());
        localBasis.evaluateFunction(ipLocal, shapeValues);

        VelocityVector v(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            // interpolate velocity at scvf
            v.axpy(shapeValues[scv.indexInElement()][0],  this->curSol()[freeFlowMomentumIdx][scv.dofIndex()]);
        }

        return v;
    }

private:
    void computeCouplingStencils_()
    {
        // TODO higher order
        const auto& momentumGridGeometry = this->problem(freeFlowMomentumIdx).gridGeometry();
        auto momentumFvGeometry = localView(momentumGridGeometry);
        this->massAndEnergyToMomentumStencils_.clear();
        this->massAndEnergyToMomentumStencils_.resize(momentumGridGeometry.gridView().size(0));

        this->momentumToMassAndEnergyStencils_.clear();
        this->momentumToMassAndEnergyStencils_.resize(momentumGridGeometry.numScv());

        for (const auto& element : elements(momentumGridGeometry.gridView()))
        {
            const auto eIdx = momentumGridGeometry.elementMapper().index(element);
            momentumFvGeometry.bindElement(element);
            for (const auto& scv : scvs(momentumFvGeometry))
            {
                this->massAndEnergyToMomentumStencils_[eIdx].push_back(scv.dofIndex());
                this->momentumToMassAndEnergyStencils_[scv.index()].push_back(eIdx);

                // extend the stencil for fluids with variable viscosity and density,
                if constexpr (FluidSystem::isCompressible(0/*phaseIdx*/))
                // if constexpr (FluidSystem::isCompressible(0/*phaseIdx*/) || !FluidSystem::viscosityIsConstant(0/*phaseIdx*/)) // TODO fix on master
                {
                    for (const auto& scvf : scvfs(momentumFvGeometry, scv))
                    {
                        if (!scvf.boundary())
                        {
                            const auto& outsideScv = momentumFvGeometry.scv(scvf.outsideScvIdx());
                            this->momentumToMassAndEnergyStencils_[scv.index()].push_back(outsideScv.elementIndex());
                        }
                    }
                }
            }
        }
    }

};

} //end namespace Dumux

#endif
