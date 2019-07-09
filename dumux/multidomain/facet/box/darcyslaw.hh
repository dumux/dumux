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
 * \copydoc Dumux::BoxFacetCouplingDarcysLaw
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FACET_COUPLING_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_FACET_COUPLING_DARCYS_LAW_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/geometry/geometryintersection.hh>
#include <dumux/common/geometry/diameter.hh>

#include <dumux/discretization/method.hh>
#include <dumux/flux/box/darcyslaw.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Darcy's law for the box scheme scheme in the context of coupled models
 *        where coupling occurs across the facets of the bulk domain elements
 *        with a lower-dimensional domain living on these facets.
 */
template<class Scalar, class FVGridGeometry>
class BoxFacetCouplingDarcysLaw
{
    using DefaultBoxDarcysLaw = BoxDarcysLaw<Scalar, FVGridGeometry>;

    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

public:

    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarCache)
    {
        // if this scvf is not on an interior boundary, use the standard law
        if (!scvf.interiorBoundary())
            return DefaultBoxDarcysLaw::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarCache);

        static const Scalar xi = getParamFromGroup<Scalar>(problem.paramGroup(), "FacetCoupling.Xi", 1.0);
        if ( !Dune::FloatCmp::eq(xi, 1.0, 1e-6) )
            DUNE_THROW(Dune::NotImplemented, "Xi != 1.0 cannot be used with the Box-Facet-Coupling scheme");

        // get some references for convenience
        const auto& fluxVarCache = elemFluxVarCache[scvf];
        const auto& shapeValues = fluxVarCache.shapeValues();
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        // Evaluate the flux in a tpfa manner using a support point inside the element.
        // Therefore, we create a segment between the integration point and another point
        // in -(K*normal) direction of the face. This segment must be long enough
        // to intersect the opposite face of the element. Then, we compute the
        // part of the segment intersecting with the element and compute the
        // point inside the element by moving half its length into the interior.
        const auto& geometry = element.geometry();
        auto distVec = mv(insideVolVars.permeability(), scvf.unitOuterNormal());
        distVec /= -1.0*distVec.two_norm();
        distVec *= diameter(geometry)*5.0; // make sure segment will be long enough

        const auto p1 = scvf.ipGlobal();
        auto p2 = scvf.ipGlobal();
        p2 += distVec;

        static constexpr int dimWorld = GridView::dimensionworld;
        using Segment = Dune::MultiLinearGeometry<Scalar, 1, dimWorld>;

        using Policy = IntersectionPolicy::SegmentPolicy<typename GridView::ctype, dimWorld>;
        using IntersectionAlgorithm = GeometryIntersection<typename Element::Geometry, Segment, Policy>;
        typename IntersectionAlgorithm::Intersection intersection;

        Segment segment(Dune::GeometryTypes::line, std::vector<GlobalPosition>({p1, p2}));
        if (!IntersectionAlgorithm::intersection(geometry, segment, intersection))
            DUNE_THROW(Dune::InvalidStateException, "Could not compute support point for flux computation");

        // use center of intersection as integration point
        auto supportPoint = intersection[0];
        supportPoint += intersection[1];
        supportPoint /= 2.0;

        std::vector< Dune::FieldVector<Scalar, 1> > supportPointShapeValues;
        fvGeometry.feLocalBasis().evaluateFunction(geometry.local(supportPoint), supportPointShapeValues);

        Scalar rho = 0.0;
        Scalar supportPressure = 0.0;
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            rho += volVars.density(phaseIdx)*shapeValues[scv.indexInElement()][0];
            supportPressure += volVars.pressure(phaseIdx)*supportPointShapeValues[scv.indexInElement()][0];
        }

        // the transmissibility on the matrix side
        const auto d = scvf.ipGlobal() - supportPoint;
        const auto tm = 1.0/d.two_norm2()*vtmv(d, insideVolVars.permeability(), scvf.unitOuterNormal());

        // compute flux depending on the user's choice of boundary types
        const auto bcTypes = problem.interiorBoundaryTypes(element, scvf);
        const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);

        // compute tpfa flux such that flux continuity holds with flux in fracture
        Scalar flux;
        if (bcTypes.hasOnlyNeumann())
        {
            // On surface grids, use sqrt of aperture as distance measur
            using std::sqrt;
            const auto df = dim == dimWorld ? 0.5*facetVolVars.extrusionFactor() : 0.5*sqrt(facetVolVars.extrusionFactor());
            const auto tf = 1.0/df*vtmv(scvf.unitOuterNormal(), facetVolVars.permeability(), scvf.unitOuterNormal());

            flux = tm*tf/(tm+tf)*(supportPressure - facetVolVars.pressure(phaseIdx))*scvf.area()*insideVolVars.extrusionFactor();
        }
        else if (bcTypes.hasOnlyDirichlet())
            flux = tm*(supportPressure - facetVolVars.pressure(phaseIdx))*scvf.area()*insideVolVars.extrusionFactor();
        // mixed boundary types are not supported
        else
            DUNE_THROW(Dune::NotImplemented, "Mixed boundary types are not supported");

        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (enableGravity)
            flux -= rho*scvf.area()*insideVolVars.extrusionFactor()
                    *vtmv(scvf.unitOuterNormal(), insideVolVars.permeability(), problem.spatialParams().gravity(scvf.center()));

        return flux;
    }

    // compute transmissibilities ti for analytical jacobians
    template<class Problem, class ElementVolumeVariables, class FluxVarCache>
    static std::vector<Scalar> calculateTransmissibilities(const Problem& problem,
                                                           const Element& element,
                                                           const FVElementGeometry& fvGeometry,
                                                           const ElementVolumeVariables& elemVolVars,
                                                           const SubControlVolumeFace& scvf,
                                                           const FluxVarCache& fluxVarCache)
    {
        DUNE_THROW(Dune::NotImplemented, "transmissibilty computation for BoxFacetCouplingDarcysLaw");
    }
};

} // end namespace Dumux

#endif
