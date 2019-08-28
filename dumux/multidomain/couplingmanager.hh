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

#ifndef DUMUX_MULTIDOMAIN_COUPLING_MANAGER_HH
#define DUMUX_MULTIDOMAIN_COUPLING_MANAGER_HH

#include <memory>
#include <tuple>
#include <vector>
#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup MultiDomain
 * \brief The interface of the coupling manager for multi domain problems
 */
template<class Traits>
class CouplingManager
{
    template<std::size_t id> using SubDomainTypeTag = typename Traits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using GridView = typename GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using ProblemWeakPtr = std::weak_ptr<const Problem<id>>;
    using Problems = typename Traits::template Tuple<ProblemWeakPtr>;

public:
    //! default type used for coupling element stencils
    template<std::size_t i, std::size_t j>
    using CouplingStencilType = std::vector<std::size_t>;

    //! the type of the solution vector
    using SolutionVector = typename Traits::SolutionVector;

    /*!
     * \name member functions concerning the coupling stencils
     */
    // \{

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
    template<std::size_t i, std::size_t j>
    const CouplingStencilType<i, j>& couplingStencil(Dune::index_constant<i> domainI,
                                                     const Element<i>& elementI,
                                                     Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "Domain i cannot be coupled to itself!");
        static_assert(AlwaysFalse<Dune::index_constant<i>>::value,
                      "The coupling manager does not implement the couplingStencil() function");
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
        curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
    }

    /*!
     * \brief update variables of domain i that depend on variables in domain j after the coupling context has been updated
     *
     * \param domainI the index of domain i
     * \param localAssemblerI the local assembler assembling the element residual of an element of domain i
     * \param elemVolVars the element volume variables (all volume variables in the element local stencil) to be updated
     * \param elemFluxVarsCache the element flux variable cache (all flux variables in the element local stencil) to be updated
     *
     * \note Such variables do not necessarily exist and then this function does nothing (default)
     * \note some examples
     *       from geomechanics: the porosity of (physical) domain i (porous medium flow) depends on the displacement vector of physical domain j (mechnanics)
     *       from domaindecomposition: the transmissibilities for fluxes of domain i to domain j depend on the permeability in domain j
     *                                 (which might depend in turn on the primary variables of domain i)
     */
    template<std::size_t i, class LocalAssemblerI, class UpdatableElementVolVars, class UpdatableFluxVarCache>
    void updateCoupledVariables(Dune::index_constant<i> domainI,
                                const LocalAssemblerI& localAssemblerI,
                                UpdatableElementVolVars& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {}

    /*!
     * \brief Updates the entire solution vector, e.g. before assembly or after grid adaption
     */
    void updateSolution(const SolutionVector& curSol)
    { curSol_ = curSol; }

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

    /*!
     * \brief set the pointers to the sub problems
     * \param problems A tuple of shared pointers to the sub problems
     */
    template<typename... SubProblems>
    void setSubProblems(const std::tuple<std::shared_ptr<SubProblems>...>& problems)
    { problems_ = problems; }

    /*!
     * \brief set a pointer to one of the sub problems
     * \param problem a pointer to the sub problem
     * \param domainIdx the domain index of the sub problem
     */
    template<class SubProblem, std::size_t i>
    void setSubProblem(std::shared_ptr<SubProblem> problem, Dune::index_constant<i> domainIdx)
    { std::get<i>(problems_) = problem; }

    /*!
     * \brief Return a reference to the sub problem
     * \param domainIdx The domain index
     */
    template<std::size_t i>
    const Problem<i>& problem(Dune::index_constant<i> domainIdx) const
    {
        if (!std::get<i>(problems_).expired())
            return *std::get<i>(problems_).lock();
        else
            DUNE_THROW(Dune::InvalidStateException, "The problem pointer was not set or has already expired. Use setSubProblems() before calling this function");
    }

protected:

    /*!
     * \brief the solution vector of the coupled problem
     * \note in case of numeric differentiation the solution vector always carries the deflected solution
     */
    SolutionVector& curSol()
    { return curSol_; }

    /*!
     * \brief the solution vector of the coupled problem
     * \note in case of numeric differentiation the solution vector always carries the deflected solution
     */
    const SolutionVector& curSol() const
    { return curSol_; }

private:
    /*!
     * \brief the solution vector of the coupled problem
     * \note in case of numeric differentiation the solution vector always carries the deflected solution
     */
    SolutionVector curSol_;

    /*!
     * \brief A tuple of std::weak_ptrs to the sub problems
     * \note these are weak pointers and not shared pointers to break the cyclic dependency between coupling manager and problems
     */
    Problems problems_;

};

} //end namespace Dumux

#endif
