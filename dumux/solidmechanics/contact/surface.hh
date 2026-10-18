// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Contact
 * \brief The boundary of a solid as a piecewise linear contact surface
 */
#ifndef DUMUX_SOLIDMECHANICS_CONTACT_SURFACE_HH
#define DUMUX_SOLIDMECHANICS_CONTACT_SURFACE_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dumux/common/math.hh>

namespace Dumux::Contact {

/*!
 * \ingroup Contact
 * \brief The boundary facets of a solid as segments (2D) or triangles (3D) between the vertex
 *        degrees of freedom
 *
 * Every facet is oriented such that its normal, \f$ (t_1, -t_0) \f$ for the segment tangent
 * \f$ \mathbf{t} = \mathbf{c}_1 - \mathbf{c}_0 \f$ in 2D and
 * \f$ (\mathbf{c}_1 - \mathbf{c}_0) \times (\mathbf{c}_2 - \mathbf{c}_0) \f$ in 3D, points out of
 * the solid in the reference configuration. Quadrilateral facets are split into two triangles.
 * Every surface node carries the reference area \f$ A_s \f$ of its dual surface patch, a share of
 * \f$ 1/d \f$ of the area of each adjacent facet.
 *
 * \tparam GridView the grid view of the solid, with dimension equal to the world dimension
 */
template<class GridView>
class ContactSurface
{
    using Scalar = typename GridView::ctype;
    static constexpr int dim = GridView::dimension;
    static_assert(dim == GridView::dimensionworld, "Contact surfaces need a solid of world dimension");
    static_assert(dim == 2 || dim == 3, "Contact surfaces are implemented in 2D and 3D");

public:
    using GlobalPosition = Dune::FieldVector<Scalar, dim>;
    //! The surface node indices of the corners of a facet
    using Facet = std::array<std::size_t, dim>;
    using Corners = std::array<GlobalPosition, dim>;
    static constexpr int dimension = dim;

    //! Unit facet normals and pseudo-normals of the surface nodes in some configuration
    struct Normals
    {
        std::vector<GlobalPosition> facet;
        std::vector<GlobalPosition> node;
    };

    //! The contact surface made of all boundary facets
    template<class DofMapper>
    ContactSurface(const GridView& gridView, const DofMapper& dofMapper)
    : ContactSurface(gridView, dofMapper, [](const auto&){ return true; })
    {}

    /*!
     * \brief The contact surface made of the boundary facets selected by a predicate
     * \param gridView the grid view of the solid
     * \param dofMapper maps vertices to the indices of the displacement degrees of freedom
     * \param isContactFacet returns whether a boundary intersection belongs to the surface
     */
    template<class DofMapper, class Selector>
    ContactSurface(const GridView& gridView, const DofMapper& dofMapper, const Selector& isContactFacet)
    {
        std::unordered_map<std::size_t, std::size_t> nodeIndex;
        const auto addNode = [&](const auto& vertex)
        {
            const auto [it, inserted] = nodeIndex.try_emplace(dofMapper.index(vertex), dofIndex_.size());
            if (inserted)
            {
                dofIndex_.push_back(it->first);
                referencePosition_.push_back(vertex.geometry().corner(0));
            }
            return it->second;
        };

        for (const auto& element : elements(gridView))
        {
            if (!element.hasBoundaryIntersections())
                continue;

            const auto refElement = Dune::referenceElement(element);
            for (const auto& intersection : intersections(gridView, element))
            {
                if (!intersection.boundary() || !isContactFacet(intersection))
                    continue;

                const int facetIdx = intersection.indexInInside();
                const int numCorners = refElement.size(facetIdx, 1, dim);
                std::array<std::size_t, 4> nodes;
                for (int i = 0; i < numCorners; ++i)
                    nodes[i] = addNode(element.template subEntity<dim>(refElement.subEntity(facetIdx, 1, i, dim)));

                const auto outerNormal = intersection.centerUnitOuterNormal();
                if (numCorners == dim)
                    addFacet_(nodes[0], nodes[1], nodes[dim-1], outerNormal);
                else if (dim == 3 && numCorners == 4)
                {
                    addFacet_(nodes[0], nodes[1], nodes[3], outerNormal);
                    addFacet_(nodes[0], nodes[3], nodes[2], outerNormal);
                }
                else
                    DUNE_THROW(Dune::NotImplemented, "Contact surface facets with " << numCorners << " corners");
            }
        }

        computeAdjacency_();
    }

    std::size_t numNodes() const
    { return dofIndex_.size(); }

    std::size_t numFacets() const
    { return facets_.size(); }

    const Facet& facet(std::size_t facetIdx) const
    { return facets_[facetIdx]; }

    //! The index of the displacement degree of freedom of a surface node
    std::size_t dofIndex(std::size_t nodeIdx) const
    { return dofIndex_[nodeIdx]; }

    const GlobalPosition& referencePosition(std::size_t nodeIdx) const
    { return referencePosition_[nodeIdx]; }

    //! The reference area of the dual surface patch of a node
    Scalar referenceArea(std::size_t nodeIdx) const
    { return referenceArea_[nodeIdx]; }

    //! The facets having the node as a corner
    const std::vector<std::size_t>& nodeFacets(std::size_t nodeIdx) const
    { return nodeFacets_[nodeIdx]; }

    //! The current positions of the surface nodes for the displacement field u
    template<class DisplacementVector>
    std::vector<GlobalPosition> positions(const DisplacementVector& u) const
    {
        std::vector<GlobalPosition> x(referencePosition_);
        for (std::size_t s = 0; s < x.size(); ++s)
            for (int i = 0; i < dim; ++i)
                x[s][i] += u[dofIndex_[s]][i];
        return x;
    }

    Corners corners(std::size_t facetIdx, const std::vector<GlobalPosition>& x) const
    {
        Corners c;
        for (int j = 0; j < dim; ++j)
            c[j] = x[facets_[facetIdx][j]];
        return c;
    }

    //! The facet normal with the length of the facet (2D) or twice its area (3D)
    static GlobalPosition normal(const Corners& c)
    {
        if constexpr (dim == 2)
            return {c[1][1] - c[0][1], c[0][0] - c[1][0]};
        else
            return crossProduct(c[1] - c[0], c[2] - c[0]);
    }

    /*!
     * \brief Facet normals and node pseudo-normals for the node positions x
     *
     * The pseudo-normal of a node is the sum of the unit normals of its facets, weighted by the
     * facet angles at the node in 3D.
     */
    Normals normals(const std::vector<GlobalPosition>& x) const
    {
        Normals n{std::vector<GlobalPosition>(numFacets()), std::vector<GlobalPosition>(numNodes(), GlobalPosition(0.0))};
        for (std::size_t f = 0; f < numFacets(); ++f)
        {
            const auto c = corners(f, x);
            n.facet[f] = normal(c);
            n.facet[f] /= n.facet[f].two_norm();
            for (int j = 0; j < dim; ++j)
            {
                if constexpr (dim == 2)
                    n.node[facets_[f][j]] += n.facet[f];
                else
                {
                    auto e1 = c[(j+1)%3] - c[j];
                    auto e2 = c[(j+2)%3] - c[j];
                    const Scalar cosAngle = std::clamp((e1*e2)/(e1.two_norm()*e2.two_norm()), Scalar(-1.0), Scalar(1.0));
                    n.node[facets_[f][j]].axpy(std::acos(cosAngle), n.facet[f]);
                }
            }
        }
        return n;
    }

private:
    void addFacet_(std::size_t n0, std::size_t n1, std::size_t n2, const GlobalPosition& outerNormal)
    {
        Facet facet;
        if constexpr (dim == 2)
            facet = {n0, n1};
        else
            facet = {n0, n1, n2};

        Corners c;
        for (int j = 0; j < dim; ++j)
            c[j] = referencePosition_[facet[j]];
        auto n = normal(c);
        if (n*outerNormal < 0.0)
        {
            std::swap(facet[dim-2], facet[dim-1]);
            n *= -1.0;
        }

        const Scalar area = dim == 2 ? n.two_norm() : 0.5*n.two_norm();
        referenceArea_.resize(referencePosition_.size(), 0.0);
        for (int j = 0; j < dim; ++j)
            referenceArea_[facet[j]] += area/dim;
        facets_.push_back(facet);
    }

    void computeAdjacency_()
    {
        nodeFacets_.assign(numNodes(), {});
        for (std::size_t f = 0; f < numFacets(); ++f)
            for (int j = 0; j < dim; ++j)
                nodeFacets_[facets_[f][j]].push_back(f);
    }

    std::vector<std::size_t> dofIndex_;
    std::vector<GlobalPosition> referencePosition_;
    std::vector<Scalar> referenceArea_;
    std::vector<Facet> facets_;
    std::vector<std::vector<std::size_t>> nodeFacets_;
};

} // end namespace Dumux::Contact

#endif
