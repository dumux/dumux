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
 * \ingroup FacetCoupling
 * \brief \copydoc Dumux::FacetCouplingPoroMechanicsCouplingManager
 */
#ifndef DUMUX_FACETCOUPLING_POROELASTIC_COUPLING_MANAGER_HH
#define DUMUX_FACETCOUPLING_POROELASTIC_COUPLING_MANAGER_HH

#include <array>
#include <algorithm>
#include <type_traits>
#include <utility>

#include <dune/common/indices.hh>
#include <dune/common/exceptions.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/method.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>
#include <dumux/geomechanics/poroelastic/couplingmanager.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Coupling manager implementation that can be used in the context
 *        of models that consider a poromechanical bulk domain (coupling
 *        between a mechanical sub-problem and a porous medium flow problem
 *        on the same grid) and a lower-dimensional porous medium flow sub-
 *        domain living on the element facets.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam BulkFacetFlowMapper Class containing maps on the coupling between
 *                             the bulk and the facet flow domain
 * \tparam BulkFacetMechMapper Class containing maps on the coupling between
 *                             the bulk mechanics and the facet flow domain
 * \tparam matrixFlowDomainId  The domain id of the bulk flow problem
 * \tparam facetFlowDomainId   The domain id of the lower-dimensional flow problem
 * \tparam mechDomainId        The domain id of the geomechanical sub-problem
 */
template< class MDTraits, class BulkFacetFlowMapper, class BulkFacetMechMapper,
          std::size_t matrixFlowDomainId = 0,
          std::size_t facetFlowDomainId = 1,
          std::size_t mechDomainId = 2,
          int lagrangeDomainId = 3>
class FacetCouplingPoroMechanicsCouplingManager
: public FacetCouplingManager< MDTraits, BulkFacetFlowMapper, matrixFlowDomainId, facetFlowDomainId >
, public FacetCouplingManager< MDTraits, BulkFacetMechMapper, mechDomainId, facetFlowDomainId >
, public PoroMechanicsCouplingManager< MDTraits, matrixFlowDomainId, mechDomainId >
{
    // convenience aliases for the underlying coupling managers
    using BulkFacetFlowManager = FacetCouplingManager< MDTraits, BulkFacetFlowMapper, matrixFlowDomainId, facetFlowDomainId >;
    using BulkFacetMechManager = FacetCouplingManager< MDTraits, BulkFacetMechMapper, mechDomainId, facetFlowDomainId >;
    using PoroMechManager = PoroMechanicsCouplingManager< MDTraits, matrixFlowDomainId, mechDomainId >;

    // domain id types
    using MatrixFlowIdType = typename MDTraits::template SubDomain<matrixFlowDomainId>::Index;
    using FacetFlowIdType = typename MDTraits::template SubDomain<facetFlowDomainId>::Index;
    using MechIdType = typename MDTraits::template SubDomain<mechDomainId>::Index;
    using LagrangeIdType = typename MDTraits::template SubDomain<lagrangeDomainId>::Index;

    // extract some types from the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Scalar = GetPropType<SubDomainTypeTag<id>, Properties::Scalar>;
    template<std::size_t id> using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using GlobalPosition = typename Element<id>::Geometry::GlobalCoordinate;
    template<std::size_t id> using GridIndexType = typename IndexTraits<GridView<id>>::GridIndex;

    // extract grid dimensions
    static constexpr int bulkDim = GridView<matrixFlowDomainId>::dimension;
    static constexpr int facetDim = GridView<facetFlowDomainId>::dimension;
    static constexpr int dimWorld = GridView<matrixFlowDomainId>::dimensionworld;

    static_assert(bulkDim == GridView<mechDomainId>::dimension,
                  "Mechanical and matrix flow domain must have same dimension");
    static_assert(bulkDim == dimWorld,
                  "Bulk dim must be equal to world dimension");
    static_assert(GridView<mechDomainId>::dimensionworld == dimWorld && GridView<facetFlowDomainId>::dimensionworld == dimWorld,
                  "World dimension must be equal for all underlying grids!");
    static_assert(GridGeometry<mechDomainId>::discMethod == DiscretizationMethod::box,
                  "This coupling manager expects the box scheme to be used in the mechanical sub-domain");
    static_assert(GridGeometry<lagrangeDomainId>::discMethod == DiscretizationMethod::fem,
                  "This coupling manager expects the lagrange domain to be discretized using finite elements");
    static_assert(std::is_same<GridView<lagrangeDomainId>, GridView<facetFlowDomainId>>::value,
                  "Facet domain and lagrange domain are expected to operate on the same grid!");

    struct ContactSurface
    {
        GridIndexType<mechDomainId> mechElemIdx; // Index of element this surface is embedded in
        GlobalPosition<mechDomainId> normal;     // Normal vector of the contact surface
        GlobalPosition<mechDomainId> tangent;    // Tangential vector of the contact surface
    };

public:
    //! export domain ids
    static constexpr auto matrixFlowId = MatrixFlowIdType();
    static constexpr auto facetFlowId = FacetFlowIdType();
    static constexpr auto mechanicsId = MechIdType();
    static constexpr auto lagrangeId = LagrangeIdType();

    //! types used for coupling stencils
    //! TODO: forward to sub-managers
    template<std::size_t i, std::size_t j>
    using CouplingStencilType = std::vector<std::size_t>;

    //! the type of the solution vector
    using SolutionVector = typename MDTraits::SolutionVector;

    //! Pull up functionalities from the parent classes
    using BulkFacetFlowManager::couplingStencil;
    using BulkFacetMechManager::couplingStencil;
    using PoroMechManager::couplingStencil;

    using BulkFacetFlowManager::isCoupled;
    using BulkFacetMechManager::isCoupled;

    using BulkFacetFlowManager::isOnInteriorBoundary;
    using BulkFacetMechManager::isOnInteriorBoundary;

    using BulkFacetFlowManager::getLowDimVolVars;
    using BulkFacetMechManager::getLowDimVolVars;

    using BulkFacetFlowManager::getLowDimElement;
    using BulkFacetMechManager::getLowDimElement;

    using BulkFacetFlowManager::getLowDimElementIndex;
    using BulkFacetMechManager::getLowDimElementIndex;

    using BulkFacetFlowManager::evalSourcesFromBulk;

    using BulkFacetFlowManager::evalCouplingResidual;
    using BulkFacetMechManager::evalCouplingResidual;
    using PoroMechManager::evalCouplingResidual;

    using BulkFacetFlowManager::extendJacobianPattern;
    using BulkFacetMechManager::extendJacobianPattern;
    using PoroMechManager::extendJacobianPattern;

    using BulkFacetFlowManager::evalAdditionalDomainDerivatives;
    using BulkFacetMechManager::evalAdditionalDomainDerivatives;
    using PoroMechManager::evalAdditionalDomainDerivatives;

    /*!
     * \brief Initialize the coupling manager.
     *
     * \param matrixFlowProblem The flow problem to be solved on the bulk domain
     * \param facetFlowProblem The flow problem to be solved on the facet domain
     * \param mechProblem The mechanical problem to be solved on the bulk domain
     * \param bulkFacetFlowMapper The mapper between the bulk and facet flow domain
     * \param bulkFacetMechMapper The mapper between the bulk mechanics and the facet flow domain
     * \tparam curSol The current solution
     */
    void init(std::shared_ptr< Problem<matrixFlowId> > matrixFlowProblem,
              std::shared_ptr< Problem<facetFlowId> > facetFlowProblem,
              std::shared_ptr< Problem<mechanicsId> > mechProblem,
              std::shared_ptr< Problem<lagrangeId> > lagrangeProblem,
              std::shared_ptr< BulkFacetFlowMapper > bulkFacetFlowMapper,
              std::shared_ptr< BulkFacetMechMapper > bulkFacetMechMapper,
              const SolutionVector& curSol)
    {
        curSol_ = curSol;
        BulkFacetFlowManager::init(matrixFlowProblem, facetFlowProblem, bulkFacetFlowMapper, curSol);
        BulkFacetMechManager::init(mechProblem, facetFlowProblem, bulkFacetMechMapper, curSol);
        PoroMechManager::init(matrixFlowProblem, mechProblem, curSol);

        // we expect the lagrange domain and facet flow domain to have identical discretizations
        if ( facetFlowProblem->fvGridGeometry().gridView().size(0)
             != lagrangeProblem->gridGeometry().gridView().size(0) )
            DUNE_THROW(Dune::InvalidStateException, "Lagrange and facet flow domain must operate on the same grid!");

        // initialize contact surfaces/stencils
        init_(bulkFacetMechMapper);
    }

    /*!
     * \brief The coupling stencil of a lagrange domain element with
     *        the facet flow domain (is empty coupling stencil).
     */
    const typename BulkFacetFlowManager::template CouplingStencilType<facetFlowId>&
    couplingStencil(LagrangeIdType domainI,
                    const Element<lagrangeId>& element,
                    FacetFlowIdType domainJ) const
    { return BulkFacetFlowManager::getEmptyStencil(domainJ); }

    /*!
     * \brief The coupling stencil of a lagrange domain element with
     *        the bulk flow domain (is empty coupling stencil).
     */
    const typename BulkFacetFlowManager::template CouplingStencilType<facetFlowId>&
    couplingStencil(LagrangeIdType domainI,
                    const Element<lagrangeId>& element,
                    MatrixFlowIdType domainJ) const
    { return BulkFacetFlowManager::getEmptyStencil(domainJ); }

    /*!
     * \brief The coupling stencil of a lagrange domain element with the mechanical domain.
     */
    const std::vector<GridIndexType<mechanicsId>>& couplingStencil(LagrangeIdType domainI,
                                                                   const Element<lagrangeId>& element,
                                                                   MechIdType domainJ) const
    {
        const auto eIdx = problem(lagrangeId).gridGeometry().elementMapper().index(element);
        return lagrangeMechCouplingStencils_[eIdx];
    }

    /*!
     * \brief The coupling stencil of a flow domain domain with
     *        the lagrange domain (is empty coupling stencil).
     */
    template<std::size_t id, std::enable_if_t<(id == matrixFlowId || id == facetFlowId), int> = 0>
    std::vector<GridIndexType<lagrangeId>> couplingStencil(Dune::index_constant<id> domainI,
                                                           const Element<id>& element,
                                                           LagrangeIdType domainJ) const
    { return {}; }

    /*!
     * \brief The coupling stencil of a mechanical domain element
     *        with the lagrange domain.
     */
    const std::vector<GridIndexType<lagrangeId>>& couplingStencil(MechIdType domainI,
                                                                  const Element<mechanicsId>& element,
                                                                  LagrangeIdType domainJ) const
    {
        const auto eIdx = problem(mechanicsId).fvGridGeometry().elementMapper().index(element);
        return mechLagrangeCouplingStencils_[eIdx];
    }

    /*!
     * \brief Computes the aperture of a sub-control volume within
     *        a given lower-dimensional element as a function of the
     *        actual mechanical deformation and the intial aperture.
     *
     * \param element The (d-1)-dimensional facet grid element
     * \param scv The (d-1)-dimensional scv for which the aperture is to be evaluated.
     * \param initialAperture The initial aperture of the scv
     */
    Scalar<facetFlowId> computeAperture(const Element<facetFlowId>& element,
                                        const typename GridGeometry<facetFlowId>::SubControlVolume& scv,
                                        Scalar<facetFlowId> initialAperture) const
    { return computeAperture(element, scv.center(), initialAperture); }

    /*!
     * \brief Computes the aperture within a lower-dimensional
     *        element at the given position as a function of the
     *        actual mechanical deformation and the intial aperture.
     *
     * \param element The (d-1)-dimensional facet grid element
     * \param globalPos The global position on the facet element.
     * \param initialAperture The initial aperture of the scv
     */
    Scalar<facetFlowId> computeAperture(const Element<facetFlowId>& element,
                                        const GlobalPosition<facetFlowId>& globalPos,
                                        Scalar<facetFlowId> initialAperture) const
    {
        const auto facetElemIdx = problem(facetFlowId).fvGridGeometry().elementMapper().index(element);
        const auto& contactSurfaces = contactSurfaces_[facetElemIdx];

        Scalar<facetFlowId> normalDispJump = 0.0;
        for (const auto& contactSurface : contactSurfaces)
        {
            const auto mechIdx = contactSurface.mechElemIdx;
            const auto mechElement = problem(mechanicsId).fvGridGeometry().element(mechIdx);
            const auto mechElemSol = elementSolution(mechElement, curSol_[mechanicsId], problem(mechanicsId).fvGridGeometry());

            const auto displacement = evalSolution(mechElement, mechElement.geometry(), mechElemSol, globalPos);
            normalDispJump -= displacement*contactSurface.normal;
        }

        return initialAperture + normalDispJump;
    }

    /*!
     * \brief Evaluates the coupling element residual of a mechanical domain element
     *        with respect to a dof in the lagrange domain (dofIdxGlobalJ). This consists
     *        of the fluxes, which is where the lagrange multiplier enters as a traction
     *        acting on the interior boundaries to the lagrange domain.
     */
    template< class MechDomainLocalAssembler >
    typename LocalResidual<mechanicsId>::ElementResidualVector
    evalCouplingResidual(MechIdType domainI,
                         const MechDomainLocalAssembler& mechDomainLocalAssembler,
                         LagrangeIdType domainJ,
                         GridIndexType<lagrangeId> dofIdxGlobalJ)
    {
        using ResidualType = typename LocalResidual<mechanicsId>::ElementResidualVector;
        ResidualType residual(mechDomainLocalAssembler.fvGeometry().numScv());
        residual = 0.0;

        const auto& localResidual = mechDomainLocalAssembler.localResidual();
        for (const auto& scvf : scvfs(mechDomainLocalAssembler.fvGeometry()))
            if (scvf.interiorBoundary())
                localResidual.evalFlux(residual,
                                       problem(mechanicsId),
                                       mechDomainLocalAssembler.element(),
                                       mechDomainLocalAssembler.fvGeometry(),
                                       mechDomainLocalAssembler.curElemVolVars(),
                                       mechDomainLocalAssembler.elemBcTypes(),
                                       mechDomainLocalAssembler.elemFluxVarsCache(),
                                       scvf);

        return residual;
    }

    /*!
     * \brief Evaluates the coupling element residual of any flow domain element
     *        with respect to a dof in the lagrange domain. This is zero (no coupling).
     */
    template< std::size_t id, class FlowDomainLocalAssembler,
              std::enable_if_t<id != mechanicsId, int> = 0 >
    typename LocalResidual<id>::ElementResidualVector
    evalCouplingResidual(Dune::index_constant<id> domainI,
                         const FlowDomainLocalAssembler& flowDomainLocalAssembler,
                         LagrangeIdType domainJ,
                         GridIndexType<lagrangeId> dofIdxGlobalJ)
    {
        using ResidualType = typename LocalResidual<id>::ElementResidualVector;
        ResidualType residual(flowDomainLocalAssembler.fvGeometry().numScv());
        residual = 0.0;
        return residual;
    }

    /*!
     * \brief Evaluates the coupling element residual of a lagrange domain element
     *        with respect to a flow domain. This is zero (no coupling).
     */
    template< class LagrangeDomainLocalAssembler, std::size_t id,
              std::enable_if_t<id != mechanicsId, int> = 0 >
    typename LocalResidual<lagrangeId>::ElementResidualVector
    evalCouplingResidual(LagrangeIdType domainI,
                         const LagrangeDomainLocalAssembler& lagrangeLocalAssembler,
                         Dune::index_constant<id> domainJ,
                         GridIndexType<lagrangeId> dofIdxGlobalJ)
    {
        using ResidualType = typename LocalResidual<lagrangeId>::ElementResidualVector;
        ResidualType residual(lagrangeLocalAssembler.feGeometry().feBasisLocalView().size());
        residual = 0.0;
        return residual;
    }

    /*!
     * \brief Evaluates the coupling element residual of a lagrange domain element
     *        with respect to the mechanical domain.
     */
    template< class LagrangeDomainLocalAssembler>
    typename LocalResidual<lagrangeId>::ElementResidualVector
    evalCouplingResidual(LagrangeIdType domainI,
                         const LagrangeDomainLocalAssembler& lagrangeLocalAssembler,
                         MechIdType domainJ,
                         GridIndexType<mechanicsId> dofIdxGlobalJ)
    {
        const auto& localResidual = lagrangeLocalAssembler.localResidual();

        if (LagrangeDomainLocalAssembler::isStandardGalerkin())
            return localResidual.eval(lagrangeLocalAssembler.element(),
                                      lagrangeLocalAssembler.feGeometry(),
                                      lagrangeLocalAssembler.curElemSol());
        else
            return localResidual.eval(lagrangeLocalAssembler.element(),
                                      lagrangeLocalAssembler.feGeometry(),
                                      lagrangeLocalAssembler.curElemSol(),
                                      lagrangeLocalAssembler.trialSpaceBasisLocalView());
    }

    /*!
     * \brief Evaluates the coupling element residual of a facet flow domain element
     *        with respect to a dof in the mechanical bulk domain (dofIdxGlobalJ). This
     *        This consists of both the source term (deformation might enter the bulk
     *        permeability and thus the transfer fluxes into the facet domain) and the
     *        fluxes on the facet domain itself, as the deformation changes the aperture
     *        and thus the permeabilities of the facet elements.
     */
    template< class FacetFlowLocalAssembler >
    typename LocalResidual<facetFlowId>::ElementResidualVector
    evalCouplingResidual(FacetFlowIdType,
                         const FacetFlowLocalAssembler& facetFlowLocalAssembler,
                         MechIdType domainJ,
                         GridIndexType<mechanicsId> dofIdxGlobalJ)
    {
        // make sure this is called for the element for which the context was set
        assert(BulkFacetFlowManager::lowDimCouplingContext().isSet);
        assert(problem(facetFlowId).fvGridGeometry().elementMapper().index(facetFlowLocalAssembler.element())
                                              == BulkFacetFlowManager::lowDimCouplingContext().elementIdx);

        // both fluxes and sources are afffected by the deformation
        const auto& localResidual = facetFlowLocalAssembler.localResidual();
        return localResidual.evalFluxAndSource(facetFlowLocalAssembler.element(),
                                               facetFlowLocalAssembler.fvGeometry(),
                                               facetFlowLocalAssembler.curElemVolVars(),
                                               facetFlowLocalAssembler.elemFluxVarsCache(),
                                               facetFlowLocalAssembler.elemBcTypes());
    }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual
     *        of an element of the lagrange domain (no extra data needed).
     */
    template<class Assembler>
    void bindCouplingContext(LagrangeIdType domainI,
                             const Element<lagrangeId>& elementI,
                             const Assembler& assembler)
    {}

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual
     *        of an element of the porous medium flow problem in the bulk domain.
     */
    template<class Assembler>
    void bindCouplingContext(MatrixFlowIdType domainI,
                             const Element<matrixFlowId>& elementI,
                             const Assembler& assembler)
    {
        BulkFacetFlowManager::bindCouplingContext(matrixFlowId, elementI, assembler);
        PoroMechManager::bindCouplingContext(matrixFlowId, elementI, assembler);
    }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual
     *        of an element of the porous medium flow problem in the facet domain.
     */
    template<class Assembler>
    void bindCouplingContext(FacetFlowIdType domainI,
                             const Element<facetFlowId>& elementI,
                             const Assembler& assembler)
    {
        BulkFacetFlowManager::bindCouplingContext(facetFlowId, elementI, assembler);
        BulkFacetMechManager::bindCouplingContext(facetFlowId, elementI, assembler);
    }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual
     *        of an element of the mechanical problem in the bulk domain.
     */
    template<class Assembler>
    void bindCouplingContext(MechIdType domainI,
                             const Element<mechanicsId>& elementI,
                             const Assembler& assembler)
    {
        PoroMechManager::bindCouplingContext(mechanicsId, elementI, assembler);
        BulkFacetMechManager::bindCouplingContext(mechanicsId, elementI, assembler);
    }

    /*!
     * \brief updates all data and variables that are necessary to evaluate the residual
     *        of an element of the lagrange domain after the solution in a sub-domain has
     *        been deflected.
     */
    template<class LocalAssemblerI, std::size_t id>
    void updateCouplingContext(LagrangeIdType domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<id> domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<id>& priVarsJ,
                               int pvIdxJ)
    {
        curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
    }

    /*!
     * \brief updates all data and variables that are necessary to evaluate the residual
     *        of an element of a sub-domain after the solution in the lagrange domain
     *        has been deflected. Since coupling to the lagrange domain only occurs via
     *        primary variables, we only deflect the solution here.
     */
    template< std::size_t id, class LocalAssemblerI >
    void updateCouplingContext(Dune::index_constant<id> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               LagrangeIdType domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<lagrangeId>& priVarsJ,
                               int pvIdxJ)
    {
        curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
    }

    /*!
     * \brief updates all data and variables that are necessary to evaluate the residual
     *        of an element of the porous medium flow problem in the bulk domain after
     *        the solution in the porous medium bulk domain has been deflected.
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(MatrixFlowIdType,
                               const LocalAssemblerI& localAssemblerI,
                               MatrixFlowIdType domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<matrixFlowId>& priVarsJ,
                               int pvIdxJ)
    {
        // always deflect the solution we return here in the public interface
        curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        BulkFacetFlowManager::updateCouplingContext(matrixFlowId, localAssemblerI, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
        PoroMechManager::updateCouplingContext(matrixFlowId, localAssemblerI, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    /*!
     * \brief updates all data and variables that are necessary to evaluate the residual
     *        of an element of the porous medium flow problem in the bulk domain after
     *        the solution in the mechanical bulk domain has been deflected.
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(MatrixFlowIdType,
                               const LocalAssemblerI& localAssemblerI,
                               MechIdType domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<mechanicsId>& priVarsJ,
                               int pvIdxJ)
    {
        // always deflect the solution we return here in the public interface
        curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        PoroMechManager::updateCouplingContext(matrixFlowId, localAssemblerI, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // The updated deformation might have an effect on the facet vol vars (i.e. extrusion)
        // deflect the solution in the flow coupling manager and rebind the context
        // note: complete re-bind might not be the most efficient solution here
        BulkFacetFlowManager::curSol()[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
        BulkFacetFlowManager::bindCouplingContext(matrixFlowId, localAssemblerI.element(), localAssemblerI.assembler());
    }

    /*!
     * \brief updates all data and variables that are necessary to evaluate the residual
     *        of an element of the porous medium flow problem in the bulk domain after
     *        the solution in the facet domain has been deflected.
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(MatrixFlowIdType,
                               const LocalAssemblerI& localAssemblerI,
                               FacetFlowIdType domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<facetFlowId>& priVarsJ,
                               int pvIdxJ)
    {
        // always deflect the solution we return here in the public interface
        curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        BulkFacetFlowManager::updateCouplingContext(matrixFlowId, localAssemblerI, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    /*!
     * \brief updates all data and variables that are necessary to evaluate the residual
     *        of an element of the facet flow problem in the facet domain after the
     *        solution in the facet domain has been deflected.
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(FacetFlowIdType,
                               const LocalAssemblerI& localAssemblerI,
                               FacetFlowIdType domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<facetFlowId>& priVarsJ,
                               int pvIdxJ)
    {
        // always deflect the solution we return here in the public interface
        curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        BulkFacetFlowManager::updateCouplingContext(facetFlowId, localAssemblerI, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
        BulkFacetMechManager::updateCouplingContext(facetFlowId, localAssemblerI, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    template<class LocalAssemblerI>
    void updateCouplingContext(FacetFlowIdType,
                               const LocalAssemblerI& localAssemblerI,
                               MechIdType domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<mechanicsId>& priVarsJ,
                               int pvIdxJ)
    {
        // always deflect the solution we return here in the public interface
        curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        BulkFacetMechManager::updateCouplingContext(facetFlowId, localAssemblerI, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // The deflected deformation might has an effect on the bulk permeabilities
        // as well. We thus simply deflect the solution in the bulk-facet flow
        // manager and rebind the context.
        // note: a complete rebind might not be the most efficient solution here
        BulkFacetFlowManager::curSol()[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];;
        BulkFacetFlowManager::bindCouplingContext(facetFlowId, localAssemblerI.element(), localAssemblerI.assembler());
    }

    /*!
     * \brief updates all data and variables that are necessary to evaluate the residual
     *        of an element of the facet flow problem in the facet domain after the
     *        solution in the bulk porous medium flow domain has been deflected.
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(FacetFlowIdType,
                               const LocalAssemblerI& localAssemblerI,
                               MatrixFlowIdType domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<matrixFlowId>& priVarsJ,
                               int pvIdxJ)
    {
        // always deflect the solution we return here in the public interface
        curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        BulkFacetFlowManager::updateCouplingContext(facetFlowId, localAssemblerI, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    /*!
     * \brief updates all data and variables that are necessary to evaluate the residual
     *        of an element of the mechanical problem in the bulk domain after the
     *        solution in the mechanical bulk domain has been deflected..
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(MechIdType,
                               const LocalAssemblerI& localAssemblerI,
                               MechIdType domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<mechanicsId>& priVarsJ,
                               int pvIdxJ)
    {
        // always deflect the solution we return here in the public interface
        curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        PoroMechManager::updateCouplingContext(mechanicsId, localAssemblerI, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
        BulkFacetMechManager::updateCouplingContext(mechanicsId, localAssemblerI, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    /*!
     * \brief updates all data and variables that are necessary to evaluate the residual
     *        of an element of the mechanical problem in the bulk domain after the
     *        solution in the porous medium flow bulk domain has been deflected.
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(MechIdType,
                               const LocalAssemblerI& localAssemblerI,
                               MatrixFlowIdType domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<matrixFlowId>& priVarsJ,
                               int pvIdxJ)
    {
        // always deflect the solution we return here in the public interface
        curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        PoroMechManager::updateCouplingContext(mechanicsId, localAssemblerI, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    /*!
     * \brief updates all data and variables that are necessary to evaluate the residual
     *        of an element of the mechanical problem in the bulk domain after the
     *        solution in the facet flow domain has been deflected.
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(MechIdType,
                               const LocalAssemblerI& localAssemblerI,
                               FacetFlowIdType domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<facetFlowId>& priVarsJ,
                               int pvIdxJ)
    {
        // always deflect the solution we return here in the public interface
        curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        BulkFacetMechManager::updateCouplingContext(mechanicsId, localAssemblerI, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    /*!
     * \brief update variables of the porous medium flow problem in the bulk domain
     *        that depend on variables in domain j after the coupling context has been updated
     */
    template<class LocalAssemblerI, class UpdatableElementVolVars, class UpdatableFluxVarCache>
    void updateCoupledVariables(MatrixFlowIdType domainI,
                                const LocalAssemblerI& localAssemblerI,
                                UpdatableElementVolVars& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        BulkFacetFlowManager::updateCoupledVariables(matrixFlowId, localAssemblerI, elemVolVars, elemFluxVarsCache);
        PoroMechManager::updateCoupledVariables(matrixFlowId, localAssemblerI, elemVolVars, elemFluxVarsCache);
    }

    /*!
     * \brief update variables of the porous medium flow problem in the facet domain
     *        that depend on variables in domain j after the coupling context has been updated
     */
    template<class LocalAssemblerI, class UpdatableElementVolVars, class UpdatableFluxVarCache>
    void updateCoupledVariables(FacetFlowIdType domainI,
                                const LocalAssemblerI& localAssemblerI,
                                UpdatableElementVolVars& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        BulkFacetFlowManager::updateCoupledVariables(facetFlowId, localAssemblerI, elemVolVars, elemFluxVarsCache);
        BulkFacetMechManager::updateCoupledVariables(facetFlowId, localAssemblerI, elemVolVars, elemFluxVarsCache);
    }

    /*!
     * \brief update variables of the geomechanical problem in the bulk domain
     *        that depend on variables in domain j after the coupling context has been updated
     */
    template<class LocalAssemblerI, class UpdatableElementVolVars, class UpdatableFluxVarCache>
    void updateCoupledVariables(MechIdType domainI,
                                const LocalAssemblerI& localAssemblerI,
                                UpdatableElementVolVars& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        PoroMechManager::updateCoupledVariables(mechanicsId, localAssemblerI, elemVolVars, elemFluxVarsCache);
        BulkFacetMechManager::updateCoupledVariables(mechanicsId, localAssemblerI, elemVolVars, elemFluxVarsCache);
    }

    //! The lagrange domain has no extended jacobian pattern
    template<class JacobianPattern>
    void extendJacobianPattern(LagrangeIdType domainI, JacobianPattern& pattern) const
    {}

    /*!
     * \brief Evaluate additional derivatives of the element residual of the lagrange domain with respect
     *        to dofs in the same domain that are not in the regular stencil (see extendJacobianPattern)
     * \note The lagrange domain has no extended jacobian pattern
     */
    template<class LagrangeLocalAssembler, class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDomainDerivatives(LagrangeIdType domainI,
                                         const LagrangeLocalAssembler& lagrangeLocalAssembler,
                                         const typename LagrangeLocalAssembler::LocalResidual::ElementResidualVector& res,
                                         JacobianMatrixDiagBlock& A,
                                         GridVariables& gridVariables)
    {}

    /*!
     * \brief updates the current solution. We have to do so in all sub-managers.
     */
    void updateSolution(const SolutionVector& sol)
    {
        curSol_ = sol;
        BulkFacetFlowManager::updateSolution(sol);
        BulkFacetMechManager::updateSolution(sol);
        PoroMechManager::updateSolution(sol);
    }

    /*!
     * \brief Returns the contact surface data structure for a sub-control
     *        volume face of the mechanical domain on an interior boundary.
     */
    const ContactSurface& getContactSurface(const Element<mechanicsId>& element,
                                            const typename GridGeometry<mechanicsId>::SubControlVolumeFace& scvf)
    {
        if (!scvf.interiorBoundary())
            DUNE_THROW(Dune::InvalidStateException, "Contact surfaces only defined for interior boundary faces");

        const auto eIdx = problem(mechanicsId).fvGridGeometry().elementMapper().index(element);
        const auto& elemLocalMap = mechContactSurfaceMap_.at(eIdx);
        const auto& idxPair = elemLocalMap.at(scvf.index());
        return contactSurfaces_[idxPair.first][idxPair.second];
    }

    //! Return a const reference to one of the flow problems
    template< std::size_t id, std::enable_if_t<(id != mechDomainId && id != lagrangeId), int> = 0 >
    const Problem<id>& problem(Dune::index_constant<id> domainId) const
    { return BulkFacetFlowManager::problem(domainId); }

    //! Return a const reference to the mechanical problem
    template< std::size_t id, std::enable_if_t<id == mechDomainId, int> = 0 >
    const Problem<id>& problem(Dune::index_constant<id> domainId) const
    { return PoroMechManager::problem(domainId); }

    //! Return a const reference to the lagrange problem
    template< std::size_t id, std::enable_if_t<(id == lagrangeId), int> = 0 >
    const Problem<id>& problem(Dune::index_constant<id> domainId) const
    { return *lagrangeProblemPtr_; }

    //! Return the current solution (all sub-managers have the entire solution)
    const SolutionVector& curSol() const
    { return curSol_; }

private:

    //! Sets up the coupling stencils and interfaces to the lagrange domain
    void init_(std::shared_ptr<BulkFacetMechMapper> bulkFacetMechMapper)
    {
        const auto& mechGG = problem(mechanicsId).fvGridGeometry();
        const auto& lagrangeGG = problem(lagrangeId).gridGeometry();

        mechLagrangeCouplingStencils_.resize(mechGG.gridView().size(0));
        lagrangeMechCouplingStencils_.resize(lagrangeGG.gridView().size(0));

        static constexpr auto bulkGridId = BulkFacetMechMapper::template gridId<bulkDim>();
        static constexpr auto lowDimGridId = BulkFacetMechMapper::template gridId<facetDim>();

        // one contact surface per coupled lagrange element
        const auto& couplingMap = bulkFacetMechMapper->couplingMap(lowDimGridId, bulkGridId);
        contactSurfaces_.resize(couplingMap.size());

        // we don't support uncoupled lagrange elements (yet?)
        if (couplingMap.size() != lagrangeGG.gridView().size(0))
            DUNE_THROW(Dune::InvalidStateException, "All lagrange domain elements are expected to be coupled");

        for (const auto& mapEntry : couplingMap)
        {
            const auto lagrangeIdx = mapEntry.first;
            const auto& embedments = mapEntry.second.embedments;
            auto& surfaces = contactSurfaces_[lagrangeIdx];
            auto& surface1 = surfaces[0];
            auto& surface2 = surfaces[1];

            if (embedments.size() == 1)
                DUNE_THROW(Dune::InvalidStateException, "Contact mechanics for boundary segment");
            if (embedments.size() != 2)
                DUNE_THROW(Dune::NotImplemented, "Contact mechanics on surface grids");

            const auto mechElementIdx1 = embedments[0].first;
            const auto mechElementIdx2 = embedments[1].first;
            const auto mechElement1 = mechGG.element(mechElementIdx1);
            const auto mechElement2 = mechGG.element(mechElementIdx2);
            const auto lagrangeElement = lagrangeGG.element(lagrangeIdx);

            auto mechFvGeometry1 = localView(mechGG);
            auto mechFvGeometry2 = localView(mechGG);
            auto lagrangeFeGeometry = localView(lagrangeGG);
            mechFvGeometry1.bindElement(mechElement1);
            mechFvGeometry2.bindElement(mechElement2);
            lagrangeFeGeometry.bind(lagrangeElement);

            // fill coupling stencils
            const auto& basisLocalView = lagrangeFeGeometry.feBasisLocalView();
            for (unsigned int localDofIdx = 0; localDofIdx < basisLocalView.size(); ++localDofIdx)
            {
                const auto dofIdx = basisLocalView.index(localDofIdx);
                mechLagrangeCouplingStencils_[mechElementIdx1].push_back(dofIdx);
                mechLagrangeCouplingStencils_[mechElementIdx2].push_back(dofIdx);
            }

            for (const auto& scv : scvs(mechFvGeometry1))
                lagrangeMechCouplingStencils_[lagrangeIdx].push_back(scv.dofIndex());
            for (const auto& scv : scvs(mechFvGeometry2))
                lagrangeMechCouplingStencils_[lagrangeIdx].push_back(scv.dofIndex());

            // set up contact surfaces
            surface1.mechElemIdx = mechElementIdx1;
            surface2.mechElemIdx = mechElementIdx2;

            // normal is the same for all scvfs on that element
            const auto& scvfList1 = embedments[0].second;
            assert(scvfList1.size() > 0);
            surface1.normal = mechFvGeometry1.scvf(scvfList1[0]).unitOuterNormal();

            // compute tangent
            if (facetDim == 2)
            {
                const auto geometry = lagrangeElement.geometry();
                surface1.tangent = geometry.corner(0) - geometry.center();
                surface1.tangent /= surface1.tangent.two_norm();
            }
            else
                DUNE_THROW(Dune::NotImplemented, "Contact mechanics in 3d");

            // simply take the negatives for second surface
            surface2.mechElemIdx = mechElementIdx2;
            surface2.normal = surface1.normal;
            surface2.tangent = surface1.tangent;
            surface2.normal *= -1.0;
            surface2.tangent *= -1.0;

            const auto& scvfList2 = embedments[1].second;

            // fill map entries
            for (auto scvfIdx : scvfList1)
                mechContactSurfaceMap_[mechElementIdx1][scvfIdx] = std::make_pair(lagrangeIdx, 0);
            for (auto scvfIdx : scvfList2)
                mechContactSurfaceMap_[mechElementIdx2][scvfIdx] = std::make_pair(lagrangeIdx, 1);
        }

        // make stencils unique
        auto makeUnique = [] (auto& v)
        {
            std::sort(v.begin(), v.end());
            v.erase( std::unique(v.begin(), v.end()), v.end() );
        };

        std::for_each(mechLagrangeCouplingStencils_.begin(), mechLagrangeCouplingStencils_.end(), makeUnique);
        std::for_each(lagrangeMechCouplingStencils_.begin(), lagrangeMechCouplingStencils_.end(), makeUnique);
    }

    // Pointer to the lagrange problem
    std::shared_ptr<Problem<lagrangeId>> lagrangeProblemPtr_;

    // The stencils between mechanical and lagrange sub-domain
    std::vector<std::vector< GridIndexType<lagrangeId> >> mechLagrangeCouplingStencils_;
    std::vector<std::vector< GridIndexType<mechanicsId> >> lagrangeMechCouplingStencils_;

    // The contact surfaces
    std::vector< std::array<ContactSurface, 2> > contactSurfaces_;

    // Allows mapping of mechanical sub-control volume faces to contact surface data
    using MechIdxType = GridIndexType<mechanicsId>;
    using ContactSurfaceIdxPair = std::pair< std::size_t, unsigned int >;
    std::unordered_map< MechIdxType, std::unordered_map<MechIdxType, ContactSurfaceIdxPair> > mechContactSurfaceMap_;

    SolutionVector curSol_;
};

} // end namespace Dumux

#endif
