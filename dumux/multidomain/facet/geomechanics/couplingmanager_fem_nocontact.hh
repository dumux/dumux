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
#ifndef DUMUX_FACETCOUPLING_POROELASTIC_FEM_COUPLING_MANAGER_HH
#define DUMUX_FACETCOUPLING_POROELASTIC_FEM_COUPLING_MANAGER_HH

#include <array>
#include <algorithm>
#include <cassert>
#include <type_traits>
#include <utility>
#include <memory>

#include <dune/common/indices.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/common/geometry/diameter.hh>
#include <dumux/common/geometry/geometryintersection.hh>
#include <dumux/common/geometry/intersectspointgeometry.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/method.hh>

#include <dumux/multidomain/glue.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>
#include <dumux/geomechanics/poroelastic/couplingmanager.hh>

namespace Dumux {

//! Element index map between the mechanical and the bulk flow sub-domains
template<std::size_t poroMechId, std::size_t matrixFlowId>
class FacetCouplingPoroMechIndexMap
{
public:
    //! default constructor
    FacetCouplingPoroMechIndexMap() = default;

    //! Constructor
    FacetCouplingPoroMechIndexMap(std::vector<std::size_t>&& mechToFlowMap,
                                  std::vector<std::size_t>&& flowToMechMap)
    : mechToFlowIndexMap_(std::move(mechToFlowMap))
    , flowToMechIndexMap_(std::move(flowToMechMap))
    {}

    //! Maps an element index of mechanical domain to index in flow domain
    std::size_t map(Dune::index_constant<poroMechId> id, std::size_t eIdx) const
    { return mechToFlowIndexMap_[eIdx]; }

    //! Maps an element index of flow domain to index in mechanical domain
    std::size_t map(Dune::index_constant<matrixFlowId> id, std::size_t eIdx) const
    { return flowToMechIndexMap_[eIdx]; }

private:
    std::vector<std::size_t> mechToFlowIndexMap_;
    std::vector<std::size_t> flowToMechIndexMap_;
};

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
 * \tparam matrixFlowDomainId  The domain id of the bulk flow problem
 * \tparam facetFlowDomainId   The domain id of the lower-dimensional flow problem
 * \tparam mechDomainId        The domain id of the geomechanical sub-problem
 */
template< class MDTraits, class BulkFacetFlowMapper,
          std::size_t matrixFlowDomainId = 0,
          std::size_t facetFlowDomainId = 1,
          std::size_t mechDomainId = 2>
class FacetCouplingPoroMechanicsCouplingManager
: public FacetCouplingManager< MDTraits, BulkFacetFlowMapper, matrixFlowDomainId, facetFlowDomainId >
, public PoroMechanicsCouplingManager< MDTraits,
                                       matrixFlowDomainId,
                                       mechDomainId,
                                       FacetCouplingPoroMechIndexMap<mechDomainId, matrixFlowDomainId> >
{
    // type of index map between the mechanical and the bulk flow sub-domain
    using BulkIndexMap = FacetCouplingPoroMechIndexMap<mechDomainId, matrixFlowDomainId>;

    // convenience aliases for the underlying coupling managers
    using BulkFacetFlowManager = FacetCouplingManager< MDTraits, BulkFacetFlowMapper, matrixFlowDomainId, facetFlowDomainId >;
    using PoroMechManager = PoroMechanicsCouplingManager< MDTraits, matrixFlowDomainId, mechDomainId, BulkIndexMap >;

    // domain id types
    using MatrixFlowIdType = typename MDTraits::template SubDomain<matrixFlowDomainId>::Index;
    using FacetFlowIdType = typename MDTraits::template SubDomain<facetFlowDomainId>::Index;
    using MechIdType = typename MDTraits::template SubDomain<mechDomainId>::Index;

    // extract some types from the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Scalar = GetPropType<SubDomainTypeTag<id>, Properties::Scalar>;
    template<std::size_t id> using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;
    template<std::size_t id> using NumEqVector = GetPropType<SubDomainTypeTag<id>, Properties::NumEqVector>;
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

    //! Volume variables in the facet flow domain
    using FacetFlowGridVariables = GetPropType<SubDomainTypeTag<facetFlowDomainId>, Properties::GridVariables>;
    using FacetFlowGridVolVars = typename FacetFlowGridVariables::GridVolumeVariables;
    using FacetFlowElemVolVars  = typename FacetFlowGridVolVars::LocalView;
    using FacetFlowVolumeVariables = typename FacetFlowGridVolVars::VolumeVariables;

public:
    //! export domain ids
    static constexpr auto matrixFlowId = MatrixFlowIdType();
    static constexpr auto facetFlowId = FacetFlowIdType();
    static constexpr auto mechanicsId = MechIdType();

    //! types used for coupling stencils
    template<std::size_t i, std::size_t j>
    using CouplingStencilType = std::vector<std::size_t>;

    //! the type of the solution vector
    using SolutionVector = typename MDTraits::SolutionVector;

    //! Pull up functionalities from the parent classes
    using BulkFacetFlowManager::couplingStencil;
    using PoroMechManager::couplingStencil;

    using BulkFacetFlowManager::evalCouplingResidual;
    using PoroMechManager::evalCouplingResidual;

    using BulkFacetFlowManager::extendJacobianPattern;
    using PoroMechManager::extendJacobianPattern;

    using BulkFacetFlowManager::evalAdditionalDomainDerivatives;
    using PoroMechManager::evalAdditionalDomainDerivatives;

    using BulkFacetFlowManager::isCoupled;
    using BulkFacetFlowManager::isOnInteriorBoundary;
    using BulkFacetFlowManager::getLowDimVolVars;
    using BulkFacetFlowManager::getLowDimElement;
    using BulkFacetFlowManager::getLowDimElementIndex;
    using BulkFacetFlowManager::evalSourcesFromBulk;

    /*!
     * \brief Initialize the coupling manager.
     *
     * \param matrixFlowProblem The flow problem to be solved on the bulk domain
     * \param facetFlowProblem The flow problem to be solved on the facet domain
     * \param mechProblem The mechanical problem to be solved on the bulk domain
     * \param bulkFacetFlowMapper The mapper between the bulk and facet flow domain
     * \tparam curSol The current solution
     */
    void init(std::shared_ptr< Problem<matrixFlowId> > matrixFlowProblem,
              std::shared_ptr< Problem<facetFlowId> > facetFlowProblem,
              std::shared_ptr< Problem<mechanicsId> > mechProblem,
              std::shared_ptr< BulkFacetFlowMapper > bulkFacetFlowMapper,
              const BulkIndexMap& bulkIndexMap,
              const SolutionVector& curSol)
    {
        bulkIndexMap_ = bulkIndexMap;
        curSol_ = curSol;

        BulkFacetFlowManager::init(matrixFlowProblem, facetFlowProblem, bulkFacetFlowMapper, curSol);
        PoroMechManager::init(matrixFlowProblem, mechProblem, curSol, bulkIndexMap);

        // initialize coupling stencils
        initMechFacetStencilsAndEmbedments_(*bulkFacetFlowMapper);
    }

    /*!
     * \brief Returns true if a mechanical domain intersection is on an interior boundary.
     */
    bool isOnInteriorBoundary(const Element<mechanicsId>& element,
                              const typename GridView<mechanicsId>::Intersection& is) const
    {
        const auto eIdx = problem(mechanicsId).gridGeometry().elementMapper().index(element);
        auto it = mechFacetEmbedmentMap_.find(eIdx);
        if (it == mechFacetEmbedmentMap_.end())
            return false;

        return std::count_if(it->second.begin(),
                             it->second.end(),
                             [&is] (const auto& p) { return p.first == is.indexInInside(); });
    }

    /*!
     * \brief Return the volume variables in the facet flow domain
     *        evaluated at a position on the boundary of the mechanical domain.
     */
    FacetFlowVolumeVariables getLowDimVolVars(const Element<mechanicsId>& element,
                                              const typename GridView<mechanicsId>::Intersection& is,
                                              const GlobalPosition<mechanicsId>& globalPos) const
    {
        const auto eIdx = problem(mechanicsId).gridGeometry().elementMapper().index(element);
        auto it1 = mechFacetEmbedmentMap_.find(eIdx);
        if (it1 == mechFacetEmbedmentMap_.end())
            DUNE_THROW(Dune::InvalidStateException, "Could not find facet flow embedment (element has no embedments)!");

        const auto& idxPairs = mechFacetEmbedmentMap_.at(eIdx);
        auto it2 = std::find_if(idxPairs.begin(),
                                idxPairs.end(),
                                [&is] (const auto& p) { return p.first == is.indexInInside(); });
        if (it2 == idxPairs.end())
            DUNE_THROW(Dune::InvalidStateException, "Could not find facet flow embedment (coinciding intersection not found)!");

        const auto facetElemIdx = it2->second;
        const auto facetElement = problem(facetFlowId).gridGeometry().element(facetElemIdx);

        auto fvGeometry = localView(problem(facetFlowId).gridGeometry());
        fvGeometry.bindElement(facetElement);

        FacetFlowVolumeVariables volVars;
        FacetCoupling::makeInterpolatedVolVars(volVars, problem(facetFlowId), curSol_[facetFlowId],
                                               fvGeometry, facetElement, facetElement.geometry(), globalPos);

        return volVars;
    }

    /*!
     * \brief The coupling stencil of a mechanical domain element with
     *        the facet flow domain.
     */
    const std::vector<GridIndexType<facetFlowId>>& couplingStencil(MechIdType domainI,
                                                                   const Element<mechanicsId>& element,
                                                                   FacetFlowIdType domainJ) const
    { return mechFacetCouplingStencils_[problem(mechanicsId).gridGeometry().elementMapper().index(element)]; }

    /*!
     * \brief The coupling stencil of a facet flow domain
     *        element with the mechanical sub-domain.
     */
    const std::vector<GridIndexType<mechanicsId>>& couplingStencil(FacetFlowIdType domainI,
                                                                   const Element<facetFlowId>& element,
                                                                   MechIdType domainJ) const
    { return facetMechCouplingStencils_[problem(facetFlowId).gridGeometry().elementMapper().index(element)]; }

    /*!
     * \brief Computes the aperture of a sub-control volume within
     *        a given lower-dimensional element as a function of the
     *        actual mechanical deformation and the intial aperture.
     *
     * \param element The (d-1)-dimensional facet grid element
     * \param scv The (d-1)-dimensional scv for which the aperture is to be evaluated.
     * \param initialAperture The initial aperture of the scv
     * \param u The displacement field to compute with
     * \todo TODO: Make this independent of scv but accept global pos in element
     */
    template<class DisplacementField>
    Scalar<facetFlowId> computeAperture(const Element<facetFlowId>& element,
                                        const typename GridGeometry<facetFlowId>::SubControlVolume& scv,
                                        Scalar<facetFlowId> initialAperture,
                                        const DisplacementField& u) const
    {
        static constexpr auto bulkGridId = BulkFacetFlowMapper::template gridId<bulkDim>();
        static constexpr auto lowDimGridId = BulkFacetFlowMapper::template gridId<facetDim>();
        const auto& couplingMap = BulkFacetFlowManager::couplingMapper().couplingMap(lowDimGridId, bulkGridId);

        const auto facetFlowElemIdx = scv.elementIndex();

        // get element index of a neighbor in the bulk flow domain
        const auto& couplingEntry = couplingMap.at(facetFlowElemIdx);
        const auto& embedments = couplingEntry.embedments;
        if (embedments.size() != 2)
            DUNE_THROW(Dune::InvalidStateException, "Invalid number of embedments found");

        // get indices of the meighboring mechanical elements
        const auto bulkMechElemIdx1 = bulkIndexMap_.map(matrixFlowId, embedments[0].first);
        const auto bulkMechElemIdx2 = bulkIndexMap_.map(matrixFlowId, embedments[1].first);
        const auto bulkMechElement1 = problem(mechanicsId).gridGeometry().element(bulkMechElemIdx1);
        const auto bulkMechElement2 = problem(mechanicsId).gridGeometry().element(bulkMechElemIdx2);

        // find intersections that coincides with this facet flow element
        const auto& idxPairs1 = mechFacetEmbedmentMap_.at(bulkMechElemIdx1);
        const auto& idxPairs2 = mechFacetEmbedmentMap_.at(bulkMechElemIdx2);

        auto it1 = std::find_if(idxPairs1.begin(), idxPairs1.end(),
                                [facetFlowElemIdx] (const auto& p) { return p.second == facetFlowElemIdx; });
        auto it2 = std::find_if(idxPairs2.begin(), idxPairs2.end(),
                                [facetFlowElemIdx] (const auto& p) { return p.second == facetFlowElemIdx; });

        if (it1 == idxPairs1.end() || it2 == idxPairs2.end())
            DUNE_THROW(Dune::InvalidStateException, "Could not find coupling intersection");

        // compute displacement jump across the scv
        const auto& mechGG = problem(mechanicsId).gridGeometry();
        const auto elemSol1 = elementSolution(bulkMechElement1, u, mechGG);
        const auto elemSol2 = elementSolution(bulkMechElement2, u, mechGG);
        const auto u1 = evalSolution(bulkMechElement1, bulkMechElement1.geometry(), mechGG, elemSol1, scv.center());
        const auto u2 = evalSolution(bulkMechElement2, bulkMechElement2.geometry(), mechGG, elemSol2, scv.center());

        // get normal of intersection
        for (const auto& is : intersections(mechGG.gridView(), bulkMechElement1))
            if (is.indexInInside() == it1->first)
                return (u2-u1)*is.centerUnitOuterNormal() + initialAperture;

        DUNE_THROW(Dune::InvalidStateException, "Could not compute aperture");
    }

    /*!
     * \brief Overload of the above function defaulting to the current displacement field.
     */
    Scalar<facetFlowId> computeAperture(const Element<facetFlowId>& element,
                                        const typename GridGeometry<facetFlowId>::SubControlVolume& scv,
                                        Scalar<facetFlowId> initialAperture) const
    { return computeAperture(element, scv, initialAperture, curSol_[mechanicsId]); }

    /*!
     * \brief Evaluates the coupling element residual of a mechanical domain element
     *        with respect to a dof in the facet flow domain (dofIdxGlobalJ).
     */
    template< class MechDomainLocalAssembler >
    typename LocalResidual<mechanicsId>::ElementResidualVector
    evalCouplingResidual(MechIdType domainI,
                         const MechDomainLocalAssembler& mechDomainLocalAssembler,
                         FacetFlowIdType domainJ,
                         GridIndexType<facetFlowId> dofIdxGlobalJ)
    {
        const auto& localResidual = mechDomainLocalAssembler.localResidual();

        if (GridGeometry<mechanicsId>::isStandardGalerkin())
            return localResidual.evalNeumannSegments(mechDomainLocalAssembler.element(),
                                                     mechDomainLocalAssembler.feGeometry(),
                                                     mechDomainLocalAssembler.curElemSol());
        else
            return localResidual.evalNeumannSegments(mechDomainLocalAssembler.element(),
                                                     mechDomainLocalAssembler.feGeometry(),
                                                     mechDomainLocalAssembler.curElemSol(),
                                                     mechDomainLocalAssembler.trialSpaceBasisLocalView());
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
        // assert(BulkFacetFlowManager::lowDimCouplingContext().isSet);
        // assert(problem(facetFlowId).fvGridGeometry().elementMapper().index(facetFlowLocalAssembler.element())
        //                                       == BulkFacetFlowManager::lowDimCouplingContext().elementIdx);

        // both fluxes and sources are afffected by the deformation
        const auto& localResidual = facetFlowLocalAssembler.localResidual();
        auto res = localResidual.evalFluxAndSource(facetFlowLocalAssembler.element(),
                                                   facetFlowLocalAssembler.fvGeometry(),
                                                   facetFlowLocalAssembler.curElemVolVars(),
                                                   facetFlowLocalAssembler.elemFluxVarsCache(),
                                                   facetFlowLocalAssembler.elemBcTypes());

        // If the residual instationary, evaluate storage
        if (!localResidual.isStationary())
            res += localResidual.evalStorage(facetFlowLocalAssembler.element(),
                                             facetFlowLocalAssembler.fvGeometry(),
                                             facetFlowLocalAssembler.prevElemVolVars(),
                                             facetFlowLocalAssembler.curElemVolVars());

        return  res;
    }

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

        // The deflected deformation might has an effect on the bulk permeabilities
        // as well. We thus simply deflect the solution in the bulk-facet flow
        // manager and rebind the context.
        // note: a complete rebind might not be the most efficient solution here
        BulkFacetFlowManager::curSol()[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
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
        // only use update of the poro-mech manager which includes all variables
        PoroMechManager::updateCoupledVariables(matrixFlowId, localAssemblerI, elemVolVars, elemFluxVarsCache);
    }

    /*!
     * \brief update variables of the porous medium flow problem in the facet domain
     *        that depend on variables in domain j after the coupling context has been updated (no caching)
     */
    template<class LocalAssemblerI, class UpdatableFluxVarCache>
    void updateCoupledVariables(FacetFlowIdType domainI,
                                const LocalAssemblerI& localAssemblerI,
                                FacetFlowElemVolVars& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        // update the element volume variables to obtain the updated apertures
        elemVolVars.bind(localAssemblerI.element(), localAssemblerI.fvGeometry(), this->curSol()[facetFlowId]);
        // maybe update transmissibilities
        BulkFacetFlowManager::updateCoupledVariables(facetFlowId, localAssemblerI, elemVolVars, elemFluxVarsCache);
    }

    /*!
     * \brief update variables of the porous medium flow problem in the facet domain
     *        that depend on variables in domain j after the coupling context has been updated (no caching)
     */
    template<class LocalAssemblerI, class UpdatableFluxVarCache>
    void updateCoupledVariables(FacetFlowIdType domainI,
                                const LocalAssemblerI& localAssemblerI,
                                FacetFlowGridVolVars& gridVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    { DUNE_THROW(Dune::NotImplemented, "Caching not yet supported for poro-elastic facet coupling models"); }

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
    }

    /*!
     * \brief updates the current solution. We have to do so in all sub-managers.
     */
    void updateSolution(const SolutionVector& sol)
    {
        curSol_ = sol;
        BulkFacetFlowManager::updateSolution(sol);
        PoroMechManager::updateSolution(sol);
    }

    //! Return a const reference to one of the flow problems
    template< std::size_t id, std::enable_if_t<(id != mechDomainId), int> = 0 >
    const Problem<id>& problem(Dune::index_constant<id> domainId) const
    { return BulkFacetFlowManager::problem(domainId); }

    //! Return a const reference to the mechanical problem
    template< std::size_t id, std::enable_if_t<id == mechDomainId, int> = 0 >
    const Problem<id>& problem(Dune::index_constant<id> domainId) const
    { return PoroMechManager::problem(domainId); }

    //! Return the current solution (all sub-managers have the entire solution)
    const SolutionVector& curSol() const
    { return curSol_; }

    //! return the index map between the mechanical and the bulk flow sub-domain
    const BulkIndexMap& bulkIndexMap() const
    { return bulkIndexMap_; }

private:
    //! Sets up the coupling stencils between mechanical and facet flow domain
    void initMechFacetStencilsAndEmbedments_(const BulkFacetFlowMapper& bulkFacetFlowMapper)
    {
        static constexpr auto bulkGridId = BulkFacetFlowMapper::template gridId<bulkDim>();
        static constexpr auto facetGridId = BulkFacetFlowMapper::template gridId<facetDim>();

        const auto& bulkMechGG = problem(mechanicsId).gridGeometry();
        const auto& facetFlowGG = problem(facetFlowId).gridGeometry();

        // set up coupling stencils mechanics domain -> facet flow domain
        facetMechCouplingStencils_.resize(facetFlowGG.gridView().size(0));
        mechFacetCouplingStencils_.resize(bulkMechGG.gridView().size(0));

        const auto& couplingMap = bulkFacetFlowMapper.couplingMap(facetGridId, bulkGridId);
        for (const auto& element : elements(facetFlowGG.gridView()))
        {
            const auto eIdx = facetFlowGG.elementMapper().index(element);
            const auto& mapEntry = couplingMap.at(eIdx);
            const auto& embedments = mapEntry.embedments;

            // get dofs of this facet element
            auto fvGeometry = localView(facetFlowGG);
            fvGeometry.bindElement(element);

            std::vector<std::size_t> facetElemDofs;
            facetElemDofs.resize(fvGeometry.numScv());
            for (const auto& scv : scvs(fvGeometry))
                facetElemDofs.push_back(scv.dofIndex());

            // fill stencils with dofs in embedments
            for (const auto& embedment : embedments)
            {
                const auto bulkFlowElemIdx = embedment.first;
                const auto bulkMechElemIdx = bulkIndexMap_.map(matrixFlowId, bulkFlowElemIdx);
                const auto bulkMechElement = bulkMechGG.element(bulkMechElemIdx);

                for (auto facetDof : facetElemDofs)
                    mechFacetCouplingStencils_[bulkMechElemIdx].push_back(facetDof);

                auto feGeometry = localView(bulkMechGG);
                feGeometry.bind(bulkMechElement);
                for (unsigned int i = 0; i < feGeometry.feBasisLocalView().size(); ++i)
                    facetMechCouplingStencils_[eIdx].push_back(feGeometry.feBasisLocalView().index(i));

                // find the intersection that overlaps with this facet element
                bool found = false;
                const auto& facetGeometry = element.geometry();
                for (const auto& is : intersections(bulkMechGG.gridView(), bulkMechElement))
                {
                    using FacetElementGeometry = typename Element<facetFlowId>::Geometry;
                    using MechIntersectionGeometry = typename GridView<mechanicsId>::Intersection::Geometry;
                    using IntersectionAlgorithm = GeometryIntersection<MechIntersectionGeometry, FacetElementGeometry>;
                    typename IntersectionAlgorithm::Intersection result;

                    if (IntersectionAlgorithm::intersection(is.geometry(), facetGeometry, result))
                    {
                        const auto idxInInside = is.indexInInside();
                        auto& idxPairs = mechFacetEmbedmentMap_[bulkMechElemIdx];
                        auto it = std::find_if(idxPairs.begin(),
                                               idxPairs.end(),
                                               [idxInInside] (const auto& p) { return p.first == idxInInside; });
                        if (it != idxPairs.end())
                            DUNE_THROW(Dune::InvalidStateException, "Found the same facet embedment twice!");

                        idxPairs.push_back(std::make_pair(idxInInside, eIdx));
                        found = true;
                        break;
                    }
                }

                if (!found)
                    DUNE_THROW(Dune::InvalidStateException, "Could not find facet domain embedment for this element");
            }
        }

        // make stencils unique
        auto makeUnique = [] (auto& v)
        {
            std::sort(v.begin(), v.end());
            v.erase( std::unique(v.begin(), v.end()), v.end() );
        };

        std::for_each(mechFacetCouplingStencils_.begin(), mechFacetCouplingStencils_.end(), makeUnique);
        std::for_each(facetMechCouplingStencils_.begin(), facetMechCouplingStencils_.end(), makeUnique);
    }

    // Index map between mechanical and bulk flow subdomains
    BulkIndexMap bulkIndexMap_;

    // Maps a mechanical domain element to pairs of facet index and corresponding facet domain element idx
    using IndexPair = std::pair<std::size_t, std::size_t>;
    using MechanicsEmbedmentMap = std::unordered_map< GridIndexType<mechanicsId>, std::vector<IndexPair> >;
    MechanicsEmbedmentMap mechFacetEmbedmentMap_;

    // Coupling stencils between mechanical and facet flow domain
    std::vector<std::vector< GridIndexType<facetFlowId> >> mechFacetCouplingStencils_;
    std::vector<std::vector< GridIndexType<mechanicsId> >> facetMechCouplingStencils_;

    // Copy of solution vector
    SolutionVector curSol_;
};

} // end namespace Dumux

#endif
