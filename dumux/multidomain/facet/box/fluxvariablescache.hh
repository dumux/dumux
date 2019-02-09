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
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dumux/discretization/box/fluxvariablescache.hh>

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
        // inside the element whose orthogonal projection is the integration point on scvf
        if (scvf.interiorBoundary())
        {
            isInteriorBoundaryCache_ = true;
            const auto& geometry = element.geometry();
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());

            const auto d1 = scvf.ipGlobal() - insideScv.dofPosition();
            const auto d2 = geometry.center() - insideScv.dofPosition();

            const auto d1Norm = d1.two_norm();
            const auto d2Norm = d2.two_norm();

            using std::tan;
            using std::acos;
            const auto angle = acos( (d1*d2)/d1Norm/d2Norm );
            const auto dm = tan(angle)*d1Norm;

            ipGlobalInside_ = scvf.unitOuterNormal();
            ipGlobalInside_ *= -1.0*dm;
            ipGlobalInside_ += scvf.ipGlobal();

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
