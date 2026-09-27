// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MortarCoupling
 * \brief The trace of a subdomain on a mortar domain, and the trace kind of a subdomain.
 *
 * A new class of subdomains is given a trace kind by a type modelling Concept::MortarTrace
 * that is constructible from the subdomain grid geometry and the mortar grid geometry, and a
 * specialization of TraceTraits selecting it. The decomposition, the model and the default
 * projectors then handle it like any other kind; how the trace cells relate to the
 * subdomain's degrees of freedom is the subdomain coupling manager's business.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_TRACE_HH
#define DUMUX_MULTIDOMAIN_MORTAR_TRACE_HH

#include <concepts>
#include <cstddef>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/concepts/mortartrace_.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/io/grid/facetgridmanager.hh>

namespace Dumux::Mortar {

#ifndef DOXYGEN
namespace Detail {

//! The measure of the part of a geometry that lies in the entities of a bounding box tree
template<typename Geometry, typename EntitySet>
typename Geometry::ctype coveredMeasure(const Geometry& geometry, const BoundingBoxTree<EntitySet>& tree)
{
    using ctype = typename Geometry::ctype;
    static constexpr int dim = Geometry::mydimension;
    static constexpr int dimWorld = Geometry::coorddimension;
    ctype measure = 0.0;
    for (const auto& intersection : intersectingEntities(geometry, tree))
        // an intersection of lower dimension, where the geometry merely touches an entity,
        // has no measure on the geometry
        if (intersection.corners().size() == dim + 1)
            measure += Dune::AffineGeometry<ctype, dim, dimWorld>(
                Dune::GeometryTypes::simplex(dim), intersection.corners()
            ).volume();
    return measure;
}

} // end namespace Detail
#endif // DOXYGEN

/*!
 * \ingroup MortarCoupling
 * \brief The trace of a bulk subdomain on a mortar domain: the boundary facets of the
 *        subdomain grid that meet the mortar, as a grid of codimension one.
 *
 * Models Concept::GridMortarTrace. The trace cells are the elements of the facet grid,
 * indexed as by the entity set; per cell, the host intersections it was extracted from are
 * available through the facet grid manager interface.
 */
template<typename HostGrid, typename FacetGrid>
class FacetTrace : public FacetGridManager<HostGrid, FacetGrid>
{
    using HostElement = typename HostGrid::template Codim<0>::Entity;
    using HostIntersection = typename HostGrid::LeafGridView::Intersection;

public:
    using HostGridView = typename HostGrid::LeafGridView;
    using GridView = typename FacetGrid::LeafGridView;
    using EntitySet = GridViewGeometricEntitySet<GridView>;
    using BoundingBoxTree = Dumux::BoundingBoxTree<EntitySet>;

    /*!
     * \brief The boundary facets of the subdomain that lie on the mortar domain.
     *
     * A facet that merely touches the mortar domain is not on it. A facet the mortar
     * domain covers only in part is refused: the mortar data is projected onto the whole
     * facet, so the subdomain grid has to resolve the boundary of the mortar domain.
     */
    template<typename GridGeometry, typename MortarGridGeometry>
        requires(std::same_as<typename GridGeometry::GridView::Grid, HostGrid>
                 and requires(const MortarGridGeometry& m) { m.boundingBoxTree(); })
    FacetTrace(const GridGeometry& gridGeometry, const MortarGridGeometry& mortarGridGeometry)
    : FacetTrace(gridGeometry.gridView(), boundaryFacetsOn_(mortarGridGeometry))
    {}

    //! The facets of the host grid view the given selector picks
    template<Concept::FacetSelector<HostElement, HostIntersection> Selector>
    FacetTrace(const HostGridView& hostGridView, const Selector& selector)
    {
        this->init(hostGridView, selector);
        entitySet_ = std::make_shared<EntitySet>(gridView());
    }

    //! The facets of the host grid the given selector picks
    template<Concept::FacetSelector<HostElement, HostIntersection> Selector>
    FacetTrace(const HostGrid& hostGrid, const Selector& selector)
    : FacetTrace(hostGrid.leafGridView(), selector)
    {}

    //! The leaf grid view of the facet grid
    GridView gridView() const
    { return this->grid().leafGridView(); }

    //! The number of trace cells
    std::size_t size() const
    { return entitySet_->size(); }

    //! The trace cells as an entity set
    std::shared_ptr<const EntitySet> entitySet() const
    { return entitySet_; }

    //! A bounding box tree over the trace cells
    const BoundingBoxTree& boundingBoxTree() const
    {
        if (!boundingBoxTree_)
            boundingBoxTree_ = std::make_unique<BoundingBoxTree>(entitySet_);
        return *boundingBoxTree_;
    }

private:
    static constexpr double coverageTolerance_ = 1e-7;

    template<typename MortarGridGeometry>
    static auto boundaryFacetsOn_(const MortarGridGeometry& mortarGridGeometry)
    {
        return [&mortarGridGeometry] (const auto&, const auto& is) {
            if (!is.boundary())
                return false;
            const auto geometry = is.geometry();
            const auto measure = geometry.volume();
            const auto covered = Detail::coveredMeasure(geometry, mortarGridGeometry.boundingBoxTree());
            if (covered < coverageTolerance_*measure)
                return false;
            if (covered < (1.0 - coverageTolerance_)*measure)
                DUNE_THROW(Dune::InvalidStateException, "The boundary facet at " << geometry.center()
                           << " is covered by the mortar domain over " << covered << " of its measure "
                           << measure << ". The subdomain grid has to resolve the boundary of the mortar domain.");
            return true;
        };
    }

    std::shared_ptr<const EntitySet> entitySet_;
    mutable std::unique_ptr<BoundingBoxTree> boundingBoxTree_;
};

/*!
 * \ingroup MortarCoupling
 * \brief The kind of trace a subdomain with the given grid geometry has on a mortar grid.
 *
 * A specialization defines the member type `type`, which models Concept::MortarTrace and is
 * constructible from the subdomain grid geometry and the mortar grid geometry, yielding the
 * part of the subdomain that meets that mortar. Bulk subdomains, whose grid has one dimension
 * more than the mortar, have a FacetTrace; subdomains of lower dimension, such as networks
 * whose boundary entities meet the mortar in points, specialize this trait for their own
 * trace kind. Pairs without a specialization define no `type`, which HasTrace detects.
 */
template<typename GridGeometry, typename MortarGrid>
struct TraceTraits {};

template<typename GridGeometry, typename MortarGrid>
    requires(int(GridGeometry::GridView::dimension) == int(MortarGrid::dimension) + 1)
struct TraceTraits<GridGeometry, MortarGrid>
{
    using type = FacetTrace<typename GridGeometry::GridView::Grid, MortarGrid>;
};

/*!
 * \ingroup MortarCoupling
 * \brief True if a trace kind is defined for the given subdomain grid geometry and mortar grid.
 */
template<typename GridGeometry, typename MortarGrid>
concept HasTrace = requires { typename TraceTraits<GridGeometry, MortarGrid>::type; }
    && Concept::MortarTrace<typename TraceTraits<GridGeometry, MortarGrid>::type>;

/*!
 * \ingroup MortarCoupling
 * \brief The trace type of a subdomain with the given grid geometry on a mortar grid.
 */
template<typename GridGeometry, typename MortarGrid>
    requires HasTrace<GridGeometry, MortarGrid>
using Trace = typename TraceTraits<GridGeometry, MortarGrid>::type;

} // end namespace Dumux::Mortar

#endif
