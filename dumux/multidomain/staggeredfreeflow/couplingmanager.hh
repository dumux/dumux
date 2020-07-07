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
class StaggeredFreeFlowCouplingManager : public CouplingManager<Traits>
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
    template<std::size_t id> using ProblemWeakPtr = std::weak_ptr<const Problem<id>>;
    using Problems = typename Traits::template Tuple<ProblemWeakPtr>;
    using Scalar = typename Traits::Scalar;
    using SolutionVector = typename Traits::SolutionVector;

    using ParentType = CouplingManager<Traits>;

    using CouplingStencilType = std::vector<std::size_t>;

    using GridVariablesTuple = typename Traits::template TupleOfSharedPtr<GridVariables>;


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

    static constexpr auto pressureIdx = VolumeVariables<freeFlowMomentumIdx>::Indices::pressureIdx;

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
        gridVariables_ = gridVariables;
        this->updateSolution(curSol);

        computeCouplingStencils_();
    }

    void init(std::shared_ptr<Problem<freeFlowMomentumIdx>> momentumProblem,
              std::shared_ptr<Problem<freeFlowMassIdx>> massProblem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol,
              const SolutionVector& prevSol)
    {
        init(momentumProblem, massProblem, gridVariables, curSol);
        prevSol_ = &prevSol;
        isTransient_ = true;
    }

    // \}


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
        bindCouplingContext(Dune::index_constant<freeFlowMomentumIdx>(), element);
        const auto& scv = (*scvs(momentumCouplingContext_[0].fvGeometry).begin());

        if constexpr (!getPropValue<SubDomainTypeTag<freeFlowMomentumIdx>, Properties::NormalizePressure>())
            return momentumCouplingContext_[0].curElemVolVars[scv].pressure();
        else
            return momentumCouplingContext_[0].curElemVolVars[scv].pressure() - this->problem(freeFlowMassIdx).initial(element)[pressureIdx];
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
        bindCouplingContext(Dune::index_constant<freeFlowMomentumIdx>(), element);
        const auto& scv = (*scvs(momentumCouplingContext_[0].fvGeometry).begin());
        return considerPreviousTimeStep ? momentumCouplingContext_[0].prevElemVolVars[scv].density()
                                        : momentumCouplingContext_[0].curElemVolVars[scv].density();
    }

    /*!
     * \brief Returns the density at a given sub control volume.
     */
    Scalar density(const Element<freeFlowMomentumIdx>& element,
                   const SubControlVolume<freeFlowMomentumIdx>& scv,
                   const bool considerPreviousTimeStep = false) const
    {
        assert(!(considerPreviousTimeStep && !isTransient_));
        bindCouplingContext(Dune::index_constant<freeFlowMomentumIdx>(), element);
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
        bindCouplingContext(Dune::index_constant<freeFlowMomentumIdx>(), element);
        const auto& scv = (*scvs(momentumCouplingContext_[0].fvGeometry).begin());
        return momentumCouplingContext_[0].curElemVolVars[scv].viscosity();
    }

     /*!
     * \brief Returns the velocity at a given sub control volume face.
     */
    Scalar faceVelocity(const Element<freeFlowMassIdx>& element,
                        const SubControlVolumeFace<freeFlowMassIdx>& scvf) const
    {
        bindCouplingContext(Dune::index_constant<freeFlowMassIdx>(), element);
        const auto& scvJ = massAndEnergyCouplingContext_[0].fvGeometry.scv(scvf.index());
        return this->curSol()[freeFlowMomentumIdx][scvJ.dofIndex()];
    }

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
    CouplingStencilType couplingStencil(Dune::index_constant<freeFlowMomentumIdx> domainI,
                                        const Element<freeFlowMomentumIdx>& elementI,
                                        Dune::index_constant<freeFlowMassIdx> domainJ) const
    {
        const auto eIdx = this->problem(freeFlowMomentumIdx).gridGeometry().elementMapper().index(elementI);
        return momentumToMassAndEnergyStencils_[eIdx];
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
        const auto eIdx = this->problem(freeFlowMomentumIdx).gridGeometry().elementMapper().index(elementI);
        if (momentumCouplingContext_.empty() || momentumCouplingContext_[0].eIdx != eIdx)
        {

            auto fvGeometry = localView(this->problem(freeFlowMassIdx).gridGeometry());
            fvGeometry.bindElement(elementI);

            auto curElemVolVars = localView(gridVars_(freeFlowMassIdx).curGridVolVars());
            curElemVolVars.bindElement(elementI, fvGeometry, this->curSol());

            auto prevElemVolVars = isTransient_ ? localView(gridVars_(freeFlowMassIdx).prevGridVolVars())
                                                : localView(gridVars_(freeFlowMassIdx).curGridVolVars());

            if (isTransient_)
                prevElemVolVars.bindElement(elementI, fvGeometry, *prevSol_);

            momentumCouplingContext_.clear();
            momentumCouplingContext_.emplace_back(MomentumCouplingContext{std::move(fvGeometry), std::move(curElemVolVars), std::move(prevElemVolVars), eIdx});
        }
    }

    void bindCouplingContext(Dune::index_constant<freeFlowMassIdx> domainI,
                             const Element<freeFlowMassIdx>& elementI) const
    {
        const auto eIdx = this->problem(freeFlowMassIdx).gridGeometry().elementMapper().index(elementI);
        if (massAndEnergyCouplingContext_.empty() || massAndEnergyCouplingContext_[0].eIdx != eIdx)
        {
            auto fvGeometry = localView(this->problem(freeFlowMomentumIdx).gridGeometry());
            fvGeometry.bindElement(elementI);

            massAndEnergyCouplingContext_.clear();
            massAndEnergyCouplingContext_.emplace_back(MassAndEnergyCouplingContext{std::move(fvGeometry), eIdx});
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
            const auto& element = localAssemblerI.element();
            const auto elemSol = elementSolution(element, this->curSol()[domainJ], problem.gridGeometry());
            const auto& fvGeometry = momentumCouplingContext_[0].fvGeometry;
            const auto& scv = (*scvs(fvGeometry).begin());

            if constexpr (ElementVolumeVariables<freeFlowMassIdx>::GridVolumeVariables::cachingEnabled)
                gridVars_(freeFlowMassIdx).curGridVolVars().volVars(scv).update(std::move(elemSol), problem, element, scv);
            else
                momentumCouplingContext_[0].curElemVolVars[scv].update(std::move(elemSol), problem, element, scv);
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
        const auto& momentumGridGeometry = this->problem(freeFlowMomentumIdx).gridGeometry();
        auto momentumFvGeometry = localView(momentumGridGeometry);
        massAndEnergyToMomentumStencils_.clear();
        massAndEnergyToMomentumStencils_.resize(momentumGridGeometry.gridView().size(0));

        const auto& massAndEnergyGridGeometry = this->problem(freeFlowMassIdx).gridGeometry();
        auto massAndEnergyFvGeometry = localView(massAndEnergyGridGeometry);
        momentumToMassAndEnergyStencils_.clear();
        momentumToMassAndEnergyStencils_.resize(massAndEnergyGridGeometry.gridView().size(0));


        for (const auto& element : elements(momentumGridGeometry.gridView()))
        {
            const auto eIdx = momentumGridGeometry.elementMapper().index(element);
            momentumFvGeometry.bindElement(element);
            for (const auto& scv : scvs(momentumFvGeometry))
                massAndEnergyToMomentumStencils_[eIdx].push_back(scv.dofIndex());

            massAndEnergyFvGeometry.bindElement(element);
            momentumToMassAndEnergyStencils_[eIdx].push_back(eIdx);
            // for (const auto& scvf : scvfs(massAndEnergyFvGeometry))
            //     if (!scvf.boundary())
            //         momentumToMassAndEnergyStencils_[eIdx].push_back(scvf.outsideScvIdx());
        }
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

} //end namespace Dumux

#endif
