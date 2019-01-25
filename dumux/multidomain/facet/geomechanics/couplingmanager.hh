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

#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>
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
          std::size_t mechDomainId = 2>
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

    // extract some types from the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Scalar = GetPropType<SubDomainTypeTag<id>, Properties::Scalar>;
    template<std::size_t id> using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using FVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id> using SubControlVolume = typename FVGridGeometry<id>::SubControlVolume;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using GridIndexType = typename IndexTraits<GridView<id>>::GridIndex;

    // extract grid dimensions
    static constexpr int bulkDim = GridView<matrixFlowDomainId>::dimension;
    static constexpr int facetDim = GridView<facetFlowDomainId>::dimension;
    static constexpr int dimWorld = GridView<matrixFlowDomainId>::dimensionworld;

    // check constraints on dimension combinations
    static_assert(bulkDim == GridView<mechDomainId>::dimension, "Mechanical and matrix flow domain must have same dimension");
    static_assert(bulkDim == dimWorld, "Bulk dim must be equal to world dimension");
    static_assert(GridView<mechDomainId>::dimensionworld == dimWorld && GridView<facetFlowDomainId>::dimensionworld == dimWorld,
                  "World dimension must be equal for all underlying grids!");

    // this coupling manager expects the box scheme to be used in the mechanical domain
    static_assert(FVGridGeometry<mechDomainId>::discMethod == DiscretizationMethod::box,
                  "This coupling manager expects the box scheme to be used in the mechanical sub-domain");

    //! Structure to store info on a mechanics-facet domain interface
    struct BulkFacetInterface
    {
        // the normal of the interface pointing "into" the facet domain
        Dune::FieldVector< Scalar<facetFlowDomainId>, dimWorld > normal;
        // the mechanical domain's dofs that lie on the facet grid element
        std::vector< GridIndexType<mechDomainId> > dofList;
    };

    // for the aperture reconstruction we store the bulk-facet interface
    // informations. There are maximum two interfaces per element.
    std::vector< Dune::ReservedVector<BulkFacetInterface, 2> > facetBulkInterfaces_;

public:
    //! export domain ids
    static constexpr auto matrixFlowId = MatrixFlowIdType();
    static constexpr auto facetFlowId = FacetFlowIdType();
    static constexpr auto mechanicsId = MechIdType();

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
    void init(std::shared_ptr< Problem<matrixFlowDomainId> > matrixFlowProblem,
              std::shared_ptr< Problem<facetFlowDomainId> > facetFlowProblem,
              std::shared_ptr< Problem<mechDomainId> > mechProblem,
              std::shared_ptr< GridVariables<matrixFlowDomainId> > matrixGridVars,
              std::shared_ptr< GridVariables<facetFlowDomainId> > facetFlowGridVars,
              std::shared_ptr< GridVariables<mechDomainId> > mechGridVars,
              std::shared_ptr< BulkFacetFlowMapper > bulkFacetFlowMapper,
              std::shared_ptr< BulkFacetMechMapper > bulkFacetMechMapper,
              const SolutionVector& curSol)
    {
        curSol_ = curSol;
        BulkFacetFlowManager::init(matrixFlowProblem, facetFlowProblem, matrixGridVars, facetFlowGridVars, bulkFacetFlowMapper, curSol);
        BulkFacetMechManager::init(mechProblem, facetFlowProblem, mechGridVars, facetFlowGridVars, bulkFacetMechMapper, curSol);
        PoroMechManager::init(matrixFlowProblem, mechProblem, matrixGridVars, mechGridVars, curSol);

        // set up the map of bulk dofs coinciding with the lower-dimensional elements
        // in order to be able to compute the deformation-dependent aperture
        const auto& bulkGG = mechProblem->fvGridGeometry();
        const auto& facetGG = facetFlowProblem->fvGridGeometry();
        facetBulkInterfaces_.resize(facetGG.gridView().size(0));

        static constexpr auto bulkGridId = BulkFacetMechMapper::template gridId<bulkDim>();
        static constexpr auto lowDimGridId = BulkFacetMechMapper::template gridId<facetDim>();

        for (const auto& mapEntry : bulkFacetMechMapper->couplingMap(lowDimGridId, bulkGridId))
        {
            const auto& embedments = mapEntry.second.embedments;
            const auto numEmbedments = embedments.size();

            if (numEmbedments > 2)
                DUNE_THROW(Dune::NotImplemented, "Coupling manager for surface grids");

            auto& interfaces = facetBulkInterfaces_[mapEntry.first];
            interfaces.resize(numEmbedments);
            for (int i = 0; i < numEmbedments; ++i)
            {
                const auto bulkElementIdx = embedments[i].first;
                const auto bulkElement = bulkGG.element(bulkElementIdx);

                auto fvGeometry = localView(bulkGG);
                fvGeometry.bindElement(bulkElement);

                // add normal (same for all coupling scvfs)
                assert(embedments[i].second.size() > 0);
                const auto& scvfList = embedments[i].second;
                interfaces[i].normal = fvGeometry.scvf(scvfList[0]).unitOuterNormal();

                // add the dofs the scvfs are connected to
                for (auto scvfIdx : scvfList)
                {
                    const auto& scvf = fvGeometry.scvf(scvfIdx);
                    const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

                    assert( std::find(interfaces[i].dofList.begin(),
                                      interfaces[i].dofList.end(),
                                      scv.dofIndex()) == interfaces[i].dofList.end() );
                    interfaces[i].dofList.push_back(scv.dofIndex());
                }
            }
        }
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
    Scalar<facetFlowDomainId> computeAperture(const Element<facetFlowDomainId>& element,
                                              const SubControlVolume<facetFlowDomainId>& scv,
                                              Scalar<facetFlowDomainId> initialAperture) const
    {
        using MechIndices = typename GetPropType<SubDomainTypeTag<mechanicsId>, Properties::ModelTraits>::Indices;
        using ScalarType = Scalar<facetFlowDomainId>;

        const auto eIdx = problem(FacetFlowIdType{}).fvGridGeometry().elementMapper().index(element);
        const auto& interfaces = facetBulkInterfaces_[eIdx];

        ScalarType deltaAperture = 0.0;
        for (const auto& interFace : interfaces)
        {
            const auto& dofList = interFace.dofList;
            assert(dofList.size() == element.subEntities(facetDim));

            Dune::FieldVector<ScalarType, bulkDim> d(0.0);
            for (auto dofIdx : dofList)
                for (int dir = 0; dir < bulkDim; ++dir)
                    d[dir] += this->curSol()[mechanicsId][dofIdx][MechIndices::u(dir)];

            deltaAperture -= (d*interFace.normal)/dofList.size();
        }

        return initialAperture + deltaAperture;
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
                         MechIdType,
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
     *        of an element of the porous medium flow problem in the bulk domain.
     */
    template<class Assembler>
    void bindCouplingContext(MatrixFlowIdType,
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
    void bindCouplingContext(FacetFlowIdType,
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
    void bindCouplingContext(MechIdType,
                             const Element<mechanicsId>& elementI,
                             const Assembler& assembler)
    {
        PoroMechManager::bindCouplingContext(mechanicsId, elementI, assembler);
        BulkFacetMechManager::bindCouplingContext(mechanicsId, elementI, assembler);
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
    void updateCoupledVariables(MatrixFlowIdType,
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
    void updateCoupledVariables(FacetFlowIdType,
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
    void updateCoupledVariables(MechIdType,
                                const LocalAssemblerI& localAssemblerI,
                                UpdatableElementVolVars& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        PoroMechManager::updateCoupledVariables(mechanicsId, localAssemblerI, elemVolVars, elemFluxVarsCache);
        BulkFacetMechManager::updateCoupledVariables(mechanicsId, localAssemblerI, elemVolVars, elemFluxVarsCache);
    }

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

    //! Return a const reference to one of the flow problems
    template< std::size_t id, std::enable_if_t<id != mechDomainId, int> = 0 >
    const Problem<id>& problem(Dune::index_constant<id> domainId) const
    { return BulkFacetFlowManager::problem(domainId); }

    //! Return a const reference to the mechanical problem
    template< std::size_t id, std::enable_if_t<id == mechDomainId, int> = 0 >
    const Problem<id>& problem(Dune::index_constant<id> domainId) const
    { return PoroMechManager::problem(domainId); }

    //! Return the current solution (all sub-managers have the entire solution)
    const SolutionVector& curSol() const
    { return curSol_; }

private:
    SolutionVector curSol_;
};

} // end namespace Dumux

#endif
