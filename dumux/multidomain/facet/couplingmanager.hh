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
 * \ingroup FacetCoupling
 * \copydoc Dumux::FacetCouplingManager
 */
#ifndef DUMUX_FACETCOUPLING_MANAGER_HH
#define DUMUX_FACETCOUPLING_MANAGER_HH

#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>

namespace Dumux {
namespace FacetCoupling{

/*!
 * \ingroup FacetCoupling
 * \brief Free function that allows the creation of a volume variables object
 *        interpolated to a given position within an element. This is the standard
 *        implementation which simply interpolates the solution to the given position
 *        and then performs a volume variables update with the interpolated solution.
 *
 * \note This assumes element-wise constant parameters for the computation of secondary
 *       variables. For heteregeneous parameter distributions a default implementation
 *       cannot be defined and an adequate overload of this function has to be provided!
 * \note For cell-centered schemes this is an unnecessary overhead because all variables
 *       are constant within the cells and a volume variables update can usually be realized
 *       more efficiently. This function is mainly to be used for the box scheme!
 */
template<class VolumeVariables, class Problem, class SolutionVector, class FVGeometry>
void makeInterpolatedVolVars(VolumeVariables& volVars,
                             const Problem& problem,
                             const SolutionVector& sol,
                             const FVGeometry& fvGeometry,
                             const typename FVGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                             const typename FVGeometry::GridGeometry::GridView::template Codim<0>::Entity::Geometry& elemGeom,
                             const typename FVGeometry::GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate& pos)
{
    // interpolate solution and set it for each entry in element solution
    auto elemSol = elementSolution(element, sol, fvGeometry.gridGeometry());
    const auto centerSol = evalSolution(element, elemGeom, fvGeometry.gridGeometry(), elemSol, pos);
    for (unsigned int i = 0; i < fvGeometry.numScv(); ++i)
        elemSol[i] = centerSol;

    // Update volume variables with the interpolated solution. Note that this standard
    // implementation only works for element-wise constant parameters as we simply use
    // the first element scv for the vol var update. For heterogeneities within the element
    // or more complex models (e.g. 2p with interface solver) a corresponding overload
    // of this function has to be provided!
    volVars.update(elemSol, problem, element, *scvs(fvGeometry).begin());
}
} // end namespace FacetCoupling

/*!
 * \ingroup FacetCoupling
 * \brief Implementation for the coupling manager between two domains of dimension d
 *        and (d-1) for models considering coupling across the bulk domain element facets.
 *        The implementations are specificto the discretization method used in the bulk
 *        domain, which is extracted automatically from the grid geometry corresponding
 *        to the provided bulk domain id. Implementations for the different methods have
 *        to be provided and included at the end of this file.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam CouplingMapper Class containing maps on the coupling between dofs of different grids
 * \tparam bulkDomainId The domain id of the bulk problem
 * \tparam lowDimDomainId The domain id of the lower-dimensional problem
 * \tparam bulkDM Discretization method used in the bulk domain
 */
template< class MDTraits,
          class CouplingMapper,
          std::size_t bulkDomainId = 0,
          std::size_t lowDimDomainId = 1,
          DiscretizationMethod bulkDM = GetPropType<typename MDTraits::template SubDomain<bulkDomainId>::TypeTag, Properties::GridGeometry>::discMethod >
class FacetCouplingManager;

/*!
 * \ingroup FacetCoupling
 * \brief Class that handles the coupling between three sub-domains in models where
 *        the coupling between the two occurs across the facets of the d- and (d-1)-
 *        dimensional domains.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam CouplingMapper Class containing maps on the coupling between dofs of different grids
 * \tparam bulkDomainId The domain id of the d-dimensional bulk problem
 * \tparam facetDomainId The domain id of the (d-1)-dimensional problem living on the bulk grid facets
 * \tparam facetDomainId The domain id of the (d-2)-dimensional problem living on the bulk grid edges
 */
template< class MDTraits,
          class CouplingMapper,
          std::size_t bulkDomainId = 0,
          std::size_t facetDomainId = 1,
          std::size_t edgeDomainId = 2 >
class FacetCouplingThreeDomainManager
: public FacetCouplingManager<MDTraits, CouplingMapper, bulkDomainId, facetDomainId>
, public FacetCouplingManager<MDTraits, CouplingMapper, facetDomainId, edgeDomainId>
{
    using BulkFacetManager = FacetCouplingManager<MDTraits, CouplingMapper, bulkDomainId, facetDomainId>;
    using FacetEdgeManager = FacetCouplingManager<MDTraits, CouplingMapper, facetDomainId, edgeDomainId>;

    // convenience aliases and instances of the domain ids
    using BulkIdType = typename MDTraits::template SubDomain<bulkDomainId>::Index;
    using FacetIdType = typename MDTraits::template SubDomain<facetDomainId>::Index;
    using EdgeIdType = typename MDTraits::template SubDomain<edgeDomainId>::Index;
    static constexpr auto bulkId = BulkIdType();
    static constexpr auto facetId = FacetIdType();
    static constexpr auto edgeId = EdgeIdType();

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    // further types specific to the sub-problems
    template<std::size_t id> using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename GridGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using GridIndexType = typename IndexTraits<GridView<id>>::GridIndex;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using ElementVolumeVariables = typename GridVariables<id>::GridVolumeVariables::LocalView;
    template<std::size_t id> using ElementFluxVariablesCache = typename GridVariables<id>::GridFluxVariablesCache::LocalView;

    // helper function to check if a domain uses mpfa
    template<std::size_t id>
    static constexpr bool usesMpfa(Dune::index_constant<id> domainId)
    { return GridGeometry<domainId>::discMethod == DiscretizationMethod::ccmpfa; }

public:
    //! types used for coupling stencils
    template<std::size_t i, std::size_t j>
    using CouplingStencilType = typename std::conditional< (j == edgeDomainId),
                                                            typename FacetEdgeManager::template CouplingStencilType<i, j>,
                                                            typename BulkFacetManager::template CouplingStencilType<i, j> >::type;

    //! the type of the solution vector
    using SolutionVector = typename MDTraits::SolutionVector;

    /*!
     * \brief Initialize the coupling manager.
     *
     * \param bulkProblem The problem to be solved on the (3d) bulk domain
     * \param facetProblem The problem to be solved on the (2d) facet domain
     * \param edgeProblem The problem to be solved on the (1d) edge domain
     * \param couplingMapper The mapper object containing the connectivity between the domains
     * \param curSol The current solution
     */
    void init(std::shared_ptr< Problem<bulkId> > bulkProblem,
              std::shared_ptr< Problem<facetId> > facetProblem,
              std::shared_ptr< Problem<edgeId> > edgeProblem,
              std::shared_ptr< CouplingMapper > couplingMapper,
              const SolutionVector& curSol)
    {
        BulkFacetManager::init(bulkProblem, facetProblem, couplingMapper, curSol);
        FacetEdgeManager::init(facetProblem, edgeProblem, couplingMapper, curSol);
    }

    //! Pull up functionalities from the parent classes
    using BulkFacetManager::couplingStencil;
    using FacetEdgeManager::couplingStencil;

    using BulkFacetManager::isCoupled;
    using FacetEdgeManager::isCoupled;

    using BulkFacetManager::isOnInteriorBoundary;
    using FacetEdgeManager::isOnInteriorBoundary;

    using BulkFacetManager::getLowDimVolVars;
    using FacetEdgeManager::getLowDimVolVars;

    using BulkFacetManager::getLowDimElement;
    using FacetEdgeManager::getLowDimElement;

    using BulkFacetManager::getLowDimElementIndex;
    using FacetEdgeManager::getLowDimElementIndex;

    using BulkFacetManager::evalSourcesFromBulk;
    using FacetEdgeManager::evalSourcesFromBulk;

    using BulkFacetManager::evalCouplingResidual;
    using FacetEdgeManager::evalCouplingResidual;

    using BulkFacetManager::bindCouplingContext;
    using FacetEdgeManager::bindCouplingContext;

    using BulkFacetManager::updateCouplingContext;
    using FacetEdgeManager::updateCouplingContext;

    using BulkFacetManager::updateCoupledVariables;
    using FacetEdgeManager::updateCoupledVariables;

    using BulkFacetManager::extendJacobianPattern;
    using FacetEdgeManager::extendJacobianPattern;

    using BulkFacetManager::evalAdditionalDomainDerivatives;
    using FacetEdgeManager::evalAdditionalDomainDerivatives;

    // extension of the jacobian pattern for the facet domain only occurs
    // within the bulk-facet coupling & for mpfa being used in the bulk domain.
    template<class JacobianPattern>
    void extendJacobianPattern(FacetIdType, JacobianPattern& pattern) const
    { BulkFacetManager::extendJacobianPattern(facetId, pattern); }

    template<class FacetLocalAssembler, class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDomainDerivatives(FacetIdType,
                                         const FacetLocalAssembler& facetLocalAssembler,
                                         const typename FacetLocalAssembler::LocalResidual::ElementResidualVector& origResiduals,
                                         JacobianMatrixDiagBlock& A,
                                         GridVariables& gridVariables)
    { BulkFacetManager::evalAdditionalDomainDerivatives(facetId, facetLocalAssembler, origResiduals, A, gridVariables); }

    /*!
     * \brief The coupling stencil of the bulk with the edge domain (empty stencil).
     */
    const CouplingStencilType<bulkId, edgeId>& couplingStencil(BulkIdType domainI,
                                                               const Element<bulkId>& element,
                                                               EdgeIdType domainJ) const
    { return FacetEdgeManager::getEmptyStencil(edgeId); }

    /*!
     * \brief The coupling stencil of the edge with the bulk domain (empty stencil).
     */
    const CouplingStencilType<edgeId, bulkId>& couplingStencil(EdgeIdType domainI,
                                                               const Element<edgeId>& element,
                                                               BulkIdType domainJ) const
    { return BulkFacetManager::getEmptyStencil(bulkId); }

    /*!
     * \brief updates the current solution. We have to overload this here
     *        to avoid ambiguity and update the solution in both managers
     */
    void updateSolution(const SolutionVector& sol)
    {
        BulkFacetManager::updateSolution(sol);
        FacetEdgeManager::updateSolution(sol);
    }

    /*!
     * \brief Interface for evaluating the coupling residual between the bulk and the edge domain.
     *        This is always zero as coupling only occurs between grids of codimension one. These
     *        overloads are provided by the two parent classes. However, we need this overload in
     *        order for the overload resolution not to fail.
     */
    template<std::size_t i,
             std::size_t j,
             class LocalAssembler,
             std::enable_if_t<((i==bulkId && j==edgeId) || ((i==edgeId && j==bulkId))), int> = 0>
    typename LocalResidual<i>::ElementResidualVector
    evalCouplingResidual(Dune::index_constant<i> domainI,
                         const LocalAssembler& localAssembler,
                         Dune::index_constant<j> domainJ,
                         GridIndexType<j> dofIdxGlobalJ)
    {
        typename LocalResidual<i>::ElementResidualVector res(1);
        res = 0.0;
        return res;
    }

    /*!
     * \brief Interface for binding the coupling context for the facet domain. In this case
     *        we have to bind both the facet -> bulk and the facet -> edge coupling context.
     */
    template< class Assembler >
    void bindCouplingContext(FacetIdType, const Element<facetId>& element, const Assembler& assembler)
    {
        BulkFacetManager::bindCouplingContext(facetId, element, assembler);
        FacetEdgeManager::bindCouplingContext(facetId, element, assembler);
    }

    /*!
     * \brief Interface for updating the coupling context of the facet domain. In this case
     *        we have to update both the facet -> bulk and the facet -> edge coupling context.
     */
    template< class FacetLocalAssembler >
    void updateCouplingContext(FacetIdType domainI,
                               const FacetLocalAssembler& facetLocalAssembler,
                               FacetIdType domainJ,
                               GridIndexType<facetId> dofIdxGlobalJ,
                               const PrimaryVariables<facetId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        BulkFacetManager::updateCouplingContext(domainI, facetLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
        FacetEdgeManager::updateCouplingContext(domainI, facetLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    /*!
     * \brief Interface for updating the coupling context between the bulk and the edge domain.
     *        We do nothing here because coupling only occurs between grids of codimension one.
     *        We have to provide this overload as the overload resolution fails otherwise.
     */
    template<std::size_t i,
             std::size_t j,
             class LocalAssembler,
             std::enable_if_t<((i==bulkId && j==edgeId) || (i==edgeId && j==bulkId)), int> = 0>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssembler& localAssembler,
                               Dune::index_constant<j> domainJ,
                               GridIndexType<j> dofIdxGlobalJ,
                               const PrimaryVariables<j>& priVarsJ,
                               unsigned int pvIdxJ)
    { /*do nothing here*/ }

    /*!
     * \brief Interface for updating the local views of the facet domain after updateCouplingContext
     *        the coupling context. In this case we have to forward the both managers as the facet
     *        domain is a part in both.
     */
    template< class FacetLocalAssembler, class UpdatableElementVolVars, class UpdatableFluxVarCache>
    void updateCoupledVariables(FacetIdType domainI,
                                const FacetLocalAssembler& facetLocalAssembler,
                                UpdatableElementVolVars& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        BulkFacetManager::updateCoupledVariables(domainI, facetLocalAssembler, elemVolVars, elemFluxVarsCache);
        FacetEdgeManager::updateCoupledVariables(domainI, facetLocalAssembler, elemVolVars, elemFluxVarsCache);
    }

    //! Return a const reference to bulk or facet problem
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == facetId), int> = 0>
    const Problem<id>& problem() const { return BulkFacetManager::template problem<id>(); }

    //! Return a reference to bulk or facet problem
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == facetId), int> = 0>
    Problem<id>& problem() { return BulkFacetManager::template problem<id>(); }

    //! Return a const reference to edge problem
    template<std::size_t id, std::enable_if_t<(id == edgeId), int> = 0>
    const Problem<id>& problem() const { return FacetEdgeManager::template problem<id>(); }

    //! Return a reference to edge problem
    template<std::size_t id, std::enable_if_t<(id == edgeId), int> = 0>
    Problem<id>& problem() { return FacetEdgeManager::template problem<id>(); }
};

} // end namespace Dumux

// Here, we have to include all available implementations
#include <dumux/multidomain/facet/box/couplingmanager.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/couplingmanager.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/couplingmanager.hh>

#endif
