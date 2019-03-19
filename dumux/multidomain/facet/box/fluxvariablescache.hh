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
 * \ingroup BoxDiscretization
 * \brief Flux variables cache class for the box scheme in the context
 *        of models considering coupling across the element facets.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FACETCOUPLING_FLUXVARIABLES_CACHE_HH
#define DUMUX_DISCRETIZATION_BOX_FACETCOUPLING_FLUXVARIABLES_CACHE_HH

#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/discretization/box/fluxvariablescache.hh>
#include <dumux/common/geometry/geometryintersection.hh>
#include <dumux/common/geometry/diameter.hh>

namespace Dumux {

/*!
 * \ingroup BoxDiscretization
 * \brief Flux variables cache class for the box scheme in the context
 *        of models considering coupling across the element facets.
 *        This class does not contain any physics-/process-dependent
 *        data. It solely stores disretization-/grid-related data.
 */
template< class Scalar, class FVGridGeometry >
class BoxFacetCouplingFluxVariablesCache
: public BoxFluxVariablesCache<Scalar, FVGridGeometry>
{
    using ParentType = BoxFluxVariablesCache<Scalar, FVGridGeometry>;

    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;

    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;

public:

    //! update the cache for an scvf
    template< class Problem, class ElementVolumeVariables >
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace& scvf)
    {
        ParentType::update(problem, element, fvGeometry, elemVolVars, scvf);

        // on interior boundaries with Neumann BCs, prepare the shape values at a point
        // inside the element whose orthogonal projection on the face is the integration point on scvf
        if (scvf.interiorBoundary())
        {
            isInteriorBoundaryCache_ = true;

            // create a segment between the integration point and another point
            // in -normal direction of the face. This segment must be long enough
            // to intersect the opposite face of the element. Then, we compute the
            // part of the segment intersecting with the element and compute the
            // point inside the element by moving half its length into the interior.
            const auto& geometry = element.geometry();
            auto distVec = scvf.unitOuterNormal();
            distVec *= -1.0;
            distVec *= diameter(geometry)*5.0; // make sure to be long enough

            auto p1 = scvf.ipGlobal();
            auto p2 = scvf.ipGlobal();
            p2 += distVec;

            static constexpr int dimWorld = GridView::dimensionworld;
            using Segment = Dune::MultiLinearGeometry<Scalar, 1, dimWorld>;

            std::vector<GlobalPosition> corners({p1, p2});
            Segment segment(Dune::GeometryTypes::line, corners);

            using Intersection = GeometryIntersection<typename Element::Geometry, Segment>;
            typename Intersection::IntersectionType intersection;

            if (!Intersection::intersection(geometry, segment, intersection))
                DUNE_THROW(Dune::InvalidStateException, "Could not compute interior integration point");

            // use center of intersection as integration point
            ipGlobalInside_ = intersection[0][0];
            ipGlobalInside_ += intersection[0][1];
            ipGlobalInside_ /= 2.0;

            fvGeometry.feLocalBasis().evaluateFunction(geometry.local(ipGlobalInside_), shapeValuesInside_);
        }
    }

    //! returns the integration point inside the element for interior boundaries
    const GlobalPosition& tpfaSupportPoint() const
    { assert(isInteriorBoundaryCache_); return ipGlobalInside_; }

    //! returns the shape values at ip inside the element for interior boundaries
    const std::vector<ShapeValue>& shapeValuesAtTpfaSupportPoint() const
    { assert(isInteriorBoundaryCache_); return shapeValuesInside_; }

private:
    bool isInteriorBoundaryCache_{false};
    GlobalPosition ipGlobalInside_;
    std::vector<ShapeValue> shapeValuesInside_;
};

} // end namespace Dumux

#endif // DUMUX_DISCRETIZATION_BOX_FACETCOUPLING_FLUXVARIABLES_CACHE_HH
