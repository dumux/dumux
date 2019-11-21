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
#include <memory>
#include <type_traits>
#include <utility>

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
    template<std::size_t id> using NumEqVector = GetPropType<SubDomainTypeTag<id>, Properties::NumEqVector>;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using GlobalPosition = typename Element<id>::Geometry::GlobalCoordinate;
    template<std::size_t id> using GridIndexType = typename IndexTraits<GridView<id>>::GridIndex;

    // Glue object between the mechanical and the lagrange-multiplier subdomain
    using MechLagrangeGlue = MultiDomainGlue<GridView<mechDomainId>, GridView<lagrangeDomainId>>;
    using MechLagrangeIntersection = typename MechLagrangeGlue::Intersection;

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

    //! Class defining a sub-set of a contact surface segment.
    class ContactSurfaceSubSegment
    {
    public:
        using Geometry = typename MechLagrangeIntersection::Geometry;

        //! Constructor. Defines the geometry and the
        //! master/slave elements in the mechanical domain.
        ContactSurfaceSubSegment(const Geometry& geometry,
                                 GridIndexType<mechDomainId> masterIdx,
                                 GridIndexType<mechDomainId> slaveIdx)
        : masterSideMechElemIdx_(masterIdx)
        , slaveSideMechElemIdx_(slaveIdx)
        , geometry_(geometry)
        {}

        //! return the master side element index
        GridIndexType<mechDomainId> getMasterSideElementIndex() const
        { return masterSideMechElemIdx_; }

        //! return the slave side element index
        GridIndexType<mechDomainId> getSlaveSideElementIndex() const
        { return slaveSideMechElemIdx_; }

        //! returns true if provided element index corresponds to master side
        bool isMasterSideElement(GridIndexType<mechDomainId> eIdx) const
        { return masterSideMechElemIdx_ == eIdx; }

        //! returns true if provided element index corresponds to master side
        bool isSlaveSideElement(GridIndexType<mechDomainId> eIdx) const
        { return slaveSideMechElemIdx_ == eIdx; }

        //! return the sub-segment geometry
        const Geometry& geometry() const
        { return geometry_; }

    private:
        GridIndexType<mechDomainId> masterSideMechElemIdx_;
        GridIndexType<mechDomainId> slaveSideMechElemIdx_;
        Geometry geometry_;
    };

    //! Class defining a segment of the contact surface.
    //! For each lagrange domain element we define one such surface,
    //! which can be further sub-divided depending on the overlap with
    //! the sub-control volume faces of the adjacend mechanical domain elements.
    class ContactSurfaceSegment
    {
    public:
        using BasisVector = GlobalPosition<mechDomainId>;
        using Basis = std::array<BasisVector, dimWorld>;

        using TractionVector = Dune::FieldVector<Scalar<lagrangeDomainId>, dimWorld>;
        using SubSegment = ContactSurfaceSubSegment;

        //! set the basis of this segment
        void setBasis(Basis&& basis)
        { basis_ = std::move(basis); }

        //! returns a basis vector of the surface
        const BasisVector& getBasisVector(unsigned int direction) const
        { return basis_[direction]; }

        //! return the number of sub-segments
        std::size_t numSubSegments() const
        { return subSegments_.size(); }

        //! return the i-th sub-segment
        const ContactSurfaceSubSegment& getSubSegment(std::size_t i) const
        { assert(i < numSubSegments()); return subSegments_[i]; }

        //! add a new sub-segment
        void addSubSegment(ContactSurfaceSubSegment&& subSegment)
        { subSegments_.emplace_back(std::move(subSegment)); }

    private:
        //! Basis vectors for the traction vector at contact
        Basis basis_;
        //! Partition of this segment into sub-segments
        std::vector<ContactSurfaceSubSegment> subSegments_;
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
        lagrangeProblemPtr_ = lagrangeProblem;
        contactIntegrationOrder_ = getParam<Scalar<lagrangeId>>("ContactProblem.ContactForceIntegrationOrder");

        BulkFacetFlowManager::init(matrixFlowProblem, facetFlowProblem, bulkFacetFlowMapper, curSol);
        BulkFacetMechManager::init(mechProblem, facetFlowProblem, bulkFacetMechMapper, curSol);
        PoroMechManager::init(matrixFlowProblem, mechProblem, curSol);

        // initialize contact surfaces/stencils
        init_();
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
     * \param u The displacement field to compute with
     * \todo TODO: Make this independent of scv but accept global pos in element
     */
    template<class DisplacementField>
    Scalar<facetFlowId> computeAperture(const Element<facetFlowId>& element,
                                        const typename GridGeometry<facetFlowId>::SubControlVolume& scv,
                                        Scalar<facetFlowId> initialAperture,
                                        const DisplacementField& u) const
    {
        static constexpr auto bulkGridId = BulkFacetMechMapper::template gridId<bulkDim>();
        static constexpr auto lowDimGridId = BulkFacetMechMapper::template gridId<facetDim>();
        const auto& couplingMap = BulkFacetMechManager::couplingMapper().couplingMap(lowDimGridId, bulkGridId);

        // get element index of a neighbor in the mechanical domain
        const auto mechElemIdx = couplingMap.at(scv.elementIndex()).embedments[0].first;

        // search in all embedded contact segments if position is contained
        const auto& map = mechContactSegmentsMap_.at(mechElemIdx);
        for (const auto& scvfToSubSegments : map)
            for (const auto& subSegmentIdxPairs : scvfToSubSegments.second)
            {
                const auto& segment = contactSurfaceSegments_[subSegmentIdxPairs.first];
                const auto& subSegment = segment.getSubSegment(subSegmentIdxPairs.second);
                if (intersectsPointGeometry(element.geometry().center(), subSegment.geometry()))
                {
                    const auto lagrangeElem = problem(lagrangeId).gridGeometry().element(subSegmentIdxPairs.first);
                    return computeAperture(lagrangeElem, element.geometry().center(), initialAperture, u);
                }
            }

        DUNE_THROW(Dune::InvalidStateException, "HUE");
    }

    /*!
     * \brief Overload of the above function defaulting to the current displacement field.
     */
    Scalar<facetFlowId> computeAperture(const Element<facetFlowId>& element,
                                        const typename GridGeometry<facetFlowId>::SubControlVolume& scv,
                                        Scalar<facetFlowId> initialAperture) const
    { return computeAperture(element, scv, initialAperture, curSol_[mechanicsId]); }

    /*!
     * \brief Computes the aperture within a lower-dimensional
     *        element at the given position as a function of the
     *        actual mechanical deformation and the intial aperture.
     *
     * \param element The (d-1)-dimensional facet grid element
     * \param globalPos The global position on the facet element.
     * \param initialAperture The initial aperture of the scv
     * \param u The displacement field to compute with
     */
    template<class DisplacementField>
    Scalar<facetFlowId> computeAperture(const Element<lagrangeId>& element,
                                        const GlobalPosition<lagrangeId>& globalPos,
                                        Scalar<lagrangeId> initialAperture,
                                        const DisplacementField& u) const
    {
        const auto deltaUN = computeNormalDisplacementJump(element, globalPos, u);
        return initialAperture - deltaUN;
    }

    /*!
     * \brief Overload of the above function defaulting to the current displacement field.
     */
    Scalar<facetFlowId> computeAperture(const Element<lagrangeId>& element,
                                        const GlobalPosition<lagrangeId>& globalPos,
                                        Scalar<lagrangeId> initialAperture) const
    { return computeAperture(element, globalPos, initialAperture, curSol_[mechanicsId]); }

    /*!
     * \brief Computes the jump in displacement within a
     *        lower-dimensional element at the given position as a
     *        function of the actual mechanical deformation.
     *
     * \param element The (d-1)-dimensional facet grid element
     * \param globalPos The global position on the facet element.
     * \param u The displacement field to compute with
     */
    template<class DisplacementField>
    GlobalPosition<mechanicsId> computeDisplacementJump(const Element<lagrangeId>& element,
                                                        const GlobalPosition<lagrangeId>& globalPos,
                                                        const DisplacementField& u) const
    {
        GlobalPosition<mechanicsId> deltaU(0.0);
        const auto& segment = getContactSurfaceSegment(element);

        // iterate over sub-segments and find the one that contains globalPos
        for (unsigned int i = 0; i < segment.numSubSegments(); ++i)
        {
            const auto& subSegment = segment.getSubSegment(i);
            const auto& subSegmentGeometry = subSegment.geometry();

            if (intersectsPointGeometry(globalPos, subSegmentGeometry))
            {
                const auto masterMechIdx = subSegment.getMasterSideElementIndex();
                const auto masterMechElement = problem(mechanicsId).fvGridGeometry().element(masterMechIdx);
                const auto masterMechElemSol = elementSolution(masterMechElement, u, problem(mechanicsId).fvGridGeometry());
                deltaU += evalSolution(masterMechElement, masterMechElement.geometry(), masterMechElemSol, globalPos);

                const auto slaveMechIdx = subSegment.getSlaveSideElementIndex();
                const auto slaveMechElement = problem(mechanicsId).fvGridGeometry().element(slaveMechIdx);
                const auto slaveMechElemSol = elementSolution(slaveMechElement, u, problem(mechanicsId).fvGridGeometry());
                deltaU -= evalSolution(slaveMechElement, slaveMechElement.geometry(), slaveMechElemSol, globalPos);
                return deltaU;
            }
        }

        DUNE_THROW(Dune::InvalidStateException, "Could not find segment which contains provided globalPos");
    }

    /*!
     * \brief Overload of the above function defaulting to the current displacement field.
     */
    GlobalPosition<mechanicsId> computeDisplacementJump(const Element<lagrangeId>& element,
                                                        const GlobalPosition<lagrangeId>& globalPos) const
    { return computeDisplacementJump(element, globalPos, curSol_[mechanicsId]); }

    /*!
     * \brief Computes the jump in tangential displacement within a
     *        lower-dimensional element at the given position as a
     *        function of the actual mechanical deformation.
     *
     * \param element The (d-1)-dimensional facet grid element
     * \param globalPos The global position on the facet element.
     * \param u The displacement field to compute with
     */
    template<class DisplacementField>
    GlobalPosition<mechanicsId> computeTangentialDisplacementJump(const Element<facetFlowId>& element,
                                                                  const GlobalPosition<facetFlowId>& globalPos,
                                                                  const DisplacementField& u) const
    {
        // compute displacement jump
        const auto deltaU = computeDisplacementJump(element, globalPos, u);

        // subtract normal part of it
        const auto& normal = getContactSurfaceSegment(element).getBasisVector(dimWorld-1);
        auto deltaUN = normal;
        deltaUN *= deltaU*normal;

        return deltaU - deltaUN;
    }

    /*!
     * \brief Overload of the above function defaulting to the current displacement field.
     */
    GlobalPosition<mechanicsId> computeTangentialDisplacementJump(const Element<facetFlowId>& element,
                                                                  const GlobalPosition<facetFlowId>& globalPos) const
    { return computeTangentialDisplacementJump(element, globalPos, curSol_[mechanicsId]); }

    /*!
     * \brief Computes the jump in normal displacement within a
     *        lower-dimensional element at the given position as a
     *        function of the actual mechanical deformation.
     *
     * \param element The (d-1)-dimensional facet grid element
     * \param globalPos The global position on the facet element.
     * \param u The displacement field to compute with
     */
    template<class DisplacementField>
    Scalar<mechanicsId> computeNormalDisplacementJump(const Element<facetFlowId>& element,
                                                      const GlobalPosition<facetFlowId>& globalPos,
                                                      const DisplacementField& u) const
    {
        // compute displacement jump
        const auto deltaU = computeDisplacementJump(element, globalPos, u);

        // evaluate the part normal to the master side
        const auto& normal = getContactSurfaceSegment(element).getBasisVector(dimWorld-1);
        return deltaU*normal;
    }

    /*!
     * \brief Overload of the above function defaulting to the current displacement field.
     */
    Scalar<mechanicsId> computeNormalDisplacementJump(const Element<facetFlowId>& element,
                                                      const GlobalPosition<facetFlowId>& globalPos) const
    { return computeNormalDisplacementJump(element, globalPos, curSol_[mechanicsId]); }

    /*!
     * \brief Returns the contact force action on a sub-control
     *        volume face of the mechanical sub-domain.
     *
     * \param element The grid element of the mechanical domain
     * \param scvf The sub-control volume face of the mechanical domain
     */
    typename ContactSurfaceSegment::TractionVector
    getContactForce(const Element<mechanicsId>& element,
                    const typename GridGeometry<mechanicsId>::SubControlVolumeFace& scvf) const
    {
        const auto eIdx = problem(mechanicsId).fvGridGeometry().elementMapper().index(element);
        const auto& map = mechContactSegmentsMap_.at(eIdx);
        const auto& mapEntry = map.at(scvf.index());

        // integrate traction over all sub-segments this scvf overlaps with
        typename ContactSurfaceSegment::TractionVector force(0.0);
        for (const auto& contactIdxPair : mapEntry)
        {
            // note that segment idx = lagrange element idx
            const auto segmentIdx = contactIdxPair.first;
            const auto subSegmentIdx = contactIdxPair.second;
            const auto lagrangeElement = problem(lagrangeId).gridGeometry().element(segmentIdx);

            const auto& segment = contactSurfaceSegments_[segmentIdx];
            const auto& subSegment = segment.getSubSegment(subSegmentIdx);
            const auto& subSegmentGeometry = subSegment.geometry();

            using QuadRules = Dune::QuadratureRules<Scalar<lagrangeId>, facetDim>;
            auto rule = QuadRules::rule(subSegmentGeometry.type(), contactIntegrationOrder_);
            for (const auto& qp : rule)
            {
                const auto& ipLocal = qp.position();
                const auto& ipGlobal = subSegmentGeometry.global(ipLocal);

                auto traction = getContactTraction(lagrangeElement, ipGlobal);
                traction *= qp.weight();
                traction *= subSegmentGeometry.integrationElement(ipLocal);
                if (subSegment.isSlaveSideElement(eIdx))
                    traction *= -1.0;
                force += traction;
            }
        }

        return force;
    }

    /*!
     * \brief Returns the lagrange multiplier (traction -> force per area)
     *        at a position on a lower-dimensional lagrange-domain element.
     *
     * \param element The grid element of the facet (lagrange) domain
     * \param pos The global position on the element
     */
    typename ContactSurfaceSegment::TractionVector
    getContactTraction(const Element<lagrangeId>& element,
                       const GlobalPosition<lagrangeId>& pos) const
    {
        const auto eg = element.geometry();
        const auto& gg = problem(lagrangeId).gridGeometry();
        const auto elemSol = elementSolution(element, curSol_[lagrangeId], gg);
        return evalSolution(element, eg, gg, elemSol, pos, true);
        // auto sigma = evalSolution(element, eg, gg, elemSol, pos, true);
        // auto f = getContactSurface(element).getBasisVector(dimWorld-1);
        // f *= sigma;
        // return f;
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

        if (GridGeometry<lagrangeId>::isStandardGalerkin())
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
    template<std::size_t i, class LocalAssemblerI, std::size_t j,
             std::enable_if_t<(i == lagrangeDomainId || j == lagrangeDomainId), int> = 0>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<j> domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<j>& priVarsJ,
                               int pvIdxJ)
    {
        curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
    }

    // /*!
    //  * \brief updates all data and variables that are necessary to evaluate the residual
    //  *        of an element of a sub-domain after the solution in the lagrange domain
    //  *        has been deflected. Since coupling to the lagrange domain only occurs via
    //  *        primary variables, we only deflect the solution here.
    //  */
    // template< std::size_t id, class LocalAssemblerI >
    // void updateCouplingContext(Dune::index_constant<id> domainI,
    //                            const LocalAssemblerI& localAssemblerI,
    //                            LagrangeIdType domainJ,
    //                            std::size_t dofIdxGlobalJ,
    //                            const PrimaryVariables<lagrangeId>& priVarsJ,
    //                            int pvIdxJ)
    // {
    //     curSol_[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
    // }

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
     * \brief Returns the contact surface segment data structure defined
     *        on a lower-dimensional element of the lagrange domain.
     */
    const ContactSurfaceSegment& getContactSurfaceSegment(const Element<lagrangeId>& element) const
    {
        const auto eIdx = problem(lagrangeId).gridGeometry().elementMapper().index(element);
        return contactSurfaceSegments_[eIdx];
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
    void init_()
    {
        const auto& mechGG = problem(mechanicsId).fvGridGeometry();
        const auto& lagrangeGG = problem(lagrangeId).gridGeometry();
        const auto glue = makeGlue(mechGG, lagrangeGG);

        mechLagrangeCouplingStencils_.resize(mechGG.gridView().size(0));
        lagrangeMechCouplingStencils_.resize(lagrangeGG.gridView().size(0));

        // one contact surface segment per lagrange domain element
        contactSurfaceSegments_.resize(lagrangeGG.gridView().size(0));

        // keep track of lagrange elements for which the basis was defined already
        std::vector<bool> hasBasis(lagrangeGG.gridView().size(0), false);

        for (const auto& is : intersections(glue))
        {
            if (is.numTargetNeighbors() != 1)
                DUNE_THROW(Dune::InvalidStateException, "Only one lagrange neighbor element expected");

            const auto lagrangeElement = is.targetEntity(0);
            const auto lagrangeElemIdx = lagrangeGG.elementMapper().index(lagrangeElement);
            auto& segment = contactSurfaceSegments_[lagrangeElemIdx];

            const auto numMechNeighbors = is.numDomainNeighbors();
            if (numMechNeighbors == 1)
                DUNE_THROW(Dune::InvalidStateException, "Contact mechanics for boundary segment");
            if (numMechNeighbors != 2)
                DUNE_THROW(Dune::NotImplemented, "Contact mechanics on surface grids");

            const auto mechElement1 = is.domainEntity(0);
            const auto mechElement2 = is.domainEntity(1);
            const auto mechElementIdx1 = mechGG.elementMapper().index(mechElement1);
            const auto mechElementIdx2 = mechGG.elementMapper().index(mechElement2);

            // (maybe) define basis for this segment
            if (!hasBasis[lagrangeElemIdx])
            {
                segment.setBasis(constructSegmentBasis(lagrangeElement, mechElement1));
                hasBasis[lagrangeElemIdx] = true;
            }

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
                lagrangeMechCouplingStencils_[lagrangeElemIdx].push_back(scv.dofIndex());
            for (const auto& scv : scvs(mechFvGeometry2))
                lagrangeMechCouplingStencils_[lagrangeElemIdx].push_back(scv.dofIndex());

            // find scvfs in the two neighbors that coincide
            const auto eps = 1e-7*diameter(lagrangeElement.geometry());
            std::vector< std::pair<unsigned int, unsigned int> > scvfPairs;
            for (const auto& scvf1 : scvfs(mechFvGeometry1))
                for (const auto& scvf2 : scvfs(mechFvGeometry2))
                    if ( (scvf1.center()-scvf2.center()).two_norm() < eps )
                        scvfPairs.emplace_back(std::make_pair(scvf1.index(), scvf2.index()));

            // create sub-segments for overlaps with scvfs
            for (const auto& pair : scvfPairs)
            {
                const auto& masterScvfGeom = mechFvGeometry1.scvf(pair.first).geometry();
                auto subGeometries = getSubIntersections_(is.geometry(), masterScvfGeom);
                for (auto& subGeometry : subGeometries)
                {
                    const auto curSubSegmentIdx = segment.numSubSegments();
                    const auto d = mechElement2.geometry().center() - mechElement1.geometry().center();
                    if (d*segment.getBasisVector(dimWorld-1) > 0.0)
                        segment.addSubSegment( ContactSurfaceSubSegment(std::move(subGeometry),
                                                                        mechElementIdx1,
                                                                        mechElementIdx2) );
                    else
                        segment.addSubSegment( ContactSurfaceSubSegment(std::move(subGeometry),
                                                                        mechElementIdx2,
                                                                        mechElementIdx1) );

                    // update maps
                    mechContactSegmentsMap_[mechElementIdx1][pair.first].push_back(std::make_pair(lagrangeElemIdx, curSubSegmentIdx));
                    mechContactSegmentsMap_[mechElementIdx2][pair.second].push_back(std::make_pair(lagrangeElemIdx, curSubSegmentIdx));
                }
            }
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

    //! Constructs the basis for the traction vector defined on a lagrange element
    typename ContactSurfaceSegment::Basis constructSegmentBasis(const Element<lagrangeId>& lgElement,
                                                                const Element<mechanicsId>& masterElement)
    {
        typename ContactSurfaceSegment::Basis basis;
        const auto& lgGeometry = lgElement.geometry();
        const auto& center = lgGeometry.center();

        if (facetDim == 1)
        {
            auto& tangent = basis[0];
            auto& normal = basis[1];

            tangent = lgGeometry.corner(0) - center;
            tangent /= tangent.two_norm();

            normal[0] = -tangent[1];
            normal[1] = tangent[0];

            // normal should point out of master side
            const auto d = center-masterElement.geometry().center();
            if (normal*d < 0.0)
            {
                tangent *= -1.0;
                normal *= -1.0;
            }
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Contact mechanics in 3d");

        return basis;
    }

    //! Computes the sub-intersections of a sub-control volume
    //! face with an element intersection. Specialization for
    //! intersectionDim == 1, where a unique result is expected.
    template< class IsGeometry, class ScvfGeometry,
              std::enable_if_t<(int(IsGeometry::mydimension) == 1), int> = 0 >
    std::vector<IsGeometry> getSubIntersections_(const IsGeometry& intersectionGeometry,
                                                 const ScvfGeometry& scvfGeometry) const
    {
        using IntersectionAlgorithm = GeometryIntersection<IsGeometry, ScvfGeometry>;

        typename IntersectionAlgorithm::Intersection result;
        if (IntersectionAlgorithm::intersection(intersectionGeometry, scvfGeometry, result))
        {
            std::vector<typename IsGeometry::GlobalCoordinate> corners;
            corners.reserve(result.size());
            for (auto& c : result) corners.push_back(c);

            return {IsGeometry(Dune::GeometryTypes::simplex(1), corners)};
        }

        return {};
    }

    //! Computes the sub-intersections of a sub-control volume
    //! face with an element intersection. Specialization for
    //! intersectionDim > 1, where the result is a triangulation.
    template< class IsGeometry, class ScvfGeometry,
              std::enable_if_t<(int(IsGeometry::mydimension) > 1), int> = 0 >
    std::vector<IsGeometry> getSubIntersections_(const IsGeometry& intersectionGeometry,
                                                 const ScvfGeometry& scvfGeometry) const
    { DUNE_THROW(Dune::NotImplemented, "Sub-Intersections in 3d world"); }

    // Pointer to the lagrange problem
    std::shared_ptr<Problem<lagrangeId>> lagrangeProblemPtr_;

    // The stencils between mechanical and lagrange sub-domain
    std::vector<std::vector< GridIndexType<lagrangeId> >> mechLagrangeCouplingStencils_;
    std::vector<std::vector< GridIndexType<mechanicsId> >> lagrangeMechCouplingStencils_;

    // The contact surfaces
    std::vector< ContactSurfaceSegment > contactSurfaceSegments_;

    // A pair of contact surface segment idx + corresponding sub-segment idx
    using SegmentSubSegmentPair = std::pair<std::size_t, std::size_t>;
    // Maps an scvf to a vector of such data pairs (one entry per overlap)
    using ScvfToContactSegmentMap = std::unordered_map< GridIndexType<mechanicsId>,
                                                        std::vector<SegmentSubSegmentPair> >;

    // Maps an element of the mechanical domain to such a map
    std::unordered_map<GridIndexType<mechanicsId>, ScvfToContactSegmentMap> mechContactSegmentsMap_;

    // Copy of solution vector
    SolutionVector curSol_;

    // Integration order for contact forces
    Scalar<lagrangeId> contactIntegrationOrder_;
};

} // end namespace Dumux

#endif
