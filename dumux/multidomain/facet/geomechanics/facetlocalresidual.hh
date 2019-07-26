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
#ifndef DUMUX_FACETCOUPLING_ELASTIC_COUPLING_MANAGER_HH
#define DUMUX_FACETCOUPLING_ELASTIC_COUPLING_MANAGER_HH

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
template< class MDTraits,
          class CouplingMapper,
          std::size_t mechDomainId = 0,
          std::size_t facetFlowDomainId = 1 >
class FacetCouplingElasticCouplingManager
: public FacetCouplingManager< MDTraits, CouplingMapper, mechDomainId, facetFlowDomainId >
{
    // convenience aliases for the underlying coupling managers
    using ParentType = FacetCouplingManager< MDTraits, CouplingMapper >;

    // domain id types
    using MechIdType = typename MDTraits::template SubDomain<mechDomainId>::Index;
    using FacetFlowIdType = typename MDTraits::template SubDomain<facetFlowDomainId>::Index;

    // extract some types from the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using FVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id> using SubControlVolume = typename FVGridGeometry<id>::SubControlVolume;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using Scalar = GetPropType<SubDomainTypeTag<id>, Properties::Scalar>;
    template<std::size_t id> using GridIndexType = typename IndexTraits<GridView<id>>::GridIndex;

    // extract grid dimensions
    static constexpr int bulkDim = GridView<mechDomainId>::dimension;
    static constexpr int facetDim = GridView<facetFlowDomainId>::dimension;
    static constexpr int dimWorld = GridView<mechDomainId>::dimensionworld;

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

    /*!
     * \brief Initialize the coupling manager.
     *
     * \param mechProblem The mechanical problem to be solved on the bulk domain
     * \param facetFlowProblem The flow problem to be solved on the facet domain
     * \param couplingMapper The mapper between the bulk and facet flow domain
     * \tparam curSol The current solution
     */
    template<class SolutionVector>
    void init(std::shared_ptr< Problem<mechDomainId> > mechProblem,
              std::shared_ptr< Problem<facetFlowDomainId> > facetFlowProblem,
              std::shared_ptr< CouplingMapper > couplingMapper,
              const SolutionVector& curSol)
    {
        ParentType::init(mechProblem, facetFlowProblem, couplingMapper, curSol);

        // set up the map of bulk dofs coinciding with the lower-dimensional elements
        // in order to be able to compute the deformation-dependent aperture
        const auto& bulkGG = mechProblem->fvGridGeometry();
        const auto& facetGG = facetFlowProblem->fvGridGeometry();
        facetBulkInterfaces_.resize(facetGG.gridView().size(0));

        static constexpr auto bulkGridId = CouplingMapper::template gridId<bulkDim>();
        static constexpr auto lowDimGridId = CouplingMapper::template gridId<facetDim>();

        for (const auto& mapEntry : couplingMapper->couplingMap(lowDimGridId, bulkGridId))
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
        using MechIndices = typename GetPropType<SubDomainTypeTag<mechDomainId>, Properties::ModelTraits>::Indices;
        using ScalarType = Scalar<facetFlowDomainId>;

        const auto eIdx = this->problem(FacetFlowIdType{}).fvGridGeometry().elementMapper().index(element);
        const auto& interfaces = facetBulkInterfaces_[eIdx];

        ScalarType deltaAperture = 0.0;
        for (const auto& interFace : interfaces)
        {
            const auto& dofList = interFace.dofList;
            assert(dofList.size() == element.subEntities(facetDim));

            Dune::FieldVector<ScalarType, bulkDim> d(0.0);
            for (auto dofIdx : dofList)
                for (int dir = 0; dir < bulkDim; ++dir)
                    d[dir] += this->curSol()[MechIdType{}][dofIdx][MechIndices::u(dir)];

            deltaAperture -= (d*interFace.normal)/dofList.size();
        }

        return initialAperture + deltaAperture;
    }
};

} // end namespace Dumux

#endif
