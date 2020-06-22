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
 * \copydoc Dumux::BoxFacetCouplingFouriersLaw
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FACET_COUPLING_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_FACET_COUPLING_FOURIERS_LAW_HH

#include <cmath>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/flux/box/fourierslaw.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Fourier's law for the box scheme in the context of coupled models
 *        where coupling occurs across the facets of the bulk domain elements
 *        with a lower-dimensional domain living on these facets.
 */
template<class TypeTag>
class BoxFacetCouplingFouriersLaw
: public FouriersLawImplementation<TypeTag, DiscretizationMethod::box>
{
    using ParentType = FouriersLawImplementation<TypeTag, DiscretizationMethod::box>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

public:

    /*!
     * \brief Compute the conductive heat flux across a sub-control volume face.
     */
    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVarsCache& elemFluxVarCache)
    {
        // if this scvf is not on an interior boundary, use the standard law
        if (!scvf.interiorBoundary())
            return ParentType::flux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarCache);

        static const Scalar xi = getParamFromGroup<Scalar>(problem.paramGroup(), "FacetCoupling.Xi", 1.0);
        if ( !Dune::FloatCmp::eq(xi, 1.0, 1e-6) )
            DUNE_THROW(Dune::NotImplemented, "Xi != 1.0 cannot be used with the Box-Facet-Coupling scheme");

        // get some references for convenience
        const auto& fluxVarCache = elemFluxVarCache[scvf];
        const auto& shapeValues = fluxVarCache.shapeValues();
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        // on interior Neumann boundaries, evaluate the flux using the facet thermal conductivity
        const auto bcTypes = problem.interiorBoundaryTypes(element, scvf);
        if (bcTypes.hasOnlyNeumann())
        {
            // compute tpfa flux from integration point to facet centerline
            const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);

            // interpolate temperature to scvf integration point
            Scalar T = 0.0;
            for (const auto& scv : scvs(fvGeometry))
                T += elemVolVars[scv].temperature()*shapeValues[scv.indexInElement()][0];

            using std::sqrt;
            // If this is a surface grid, use the square root of the facet extrusion factor
            // as an approximate average distance from scvf ip to facet center
            using std::sqrt;
            const auto a = facetVolVars.extrusionFactor();
            auto gradT = scvf.unitOuterNormal();
            gradT *= dim == dimWorld ? 0.5*a : 0.5*sqrt(a);
            gradT /= gradT.two_norm2();
            gradT *= (facetVolVars.temperature() - T);

            return -1.0*Extrusion::area(scvf)
                       *insideVolVars.extrusionFactor()
                       *vtmv(scvf.unitOuterNormal(),
                             facetVolVars.effectiveThermalConductivity(),
                             gradT);
        }

        // on interior Dirichlet boundaries use the facet temperature and evaluate flux
        else if (bcTypes.hasOnlyDirichlet())
        {
            // create vector with nodal temperatures
            std::vector<Scalar> temperatures(element.subEntities(dim));
            for (const auto& scv : scvs(fvGeometry))
                temperatures[scv.localDofIndex()] = elemVolVars[scv].temperature();

            // substitute with facet temperatures for those scvs touching this facet
            for (const auto& scvfJ : scvfs(fvGeometry))
                if (scvfJ.interiorBoundary() && scvfJ.facetIndexInElement() == scvf.facetIndexInElement())
                    temperatures[ fvGeometry.scv(scvfJ.insideScvIdx()).localDofIndex() ]
                             = problem.couplingManager().getLowDimVolVars(element, scvfJ).temperature();

            // evaluate gradT at integration point
            Dune::FieldVector<Scalar, dimWorld> gradT(0.0);
            for (const auto& scv : scvs(fvGeometry))
                gradT.axpy(temperatures[scv.localDofIndex()], fluxVarCache.gradN(scv.indexInElement()));

            // apply matrix thermal conductivity and return the flux
            return -1.0*Extrusion::area(scvf)
                       *insideVolVars.extrusionFactor()
                       *vtmv(scvf.unitOuterNormal(),
                         insideVolVars.effectiveThermalConductivity(),
                         gradT);
        }

        // mixed boundary types are not supported
        else
            DUNE_THROW(Dune::NotImplemented, "Mixed boundary types are not supported");
    }

    // compute transmissibilities ti for analytical jacobians
    template<class Problem, class ElementVolumeVariables, class FluxVarCache>
    static std::vector<Scalar> calculateTransmissibilities(const Problem& problem,
                                                           const Element& element,
                                                           const FVElementGeometry& fvGeometry,
                                                           const ElementVolumeVariables& elemVolVars,
                                                           const SubControlVolumeFace& scvf,
                                                           const FluxVarCache& fluxVarCache,
                                                           unsigned int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "transmissibilty computation for BoxFacetCouplingFouriersLaw");
    }
};

} // end namespace Dumux

#endif
