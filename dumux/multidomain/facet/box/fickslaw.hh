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
 * \copydoc Dumux::BoxFacetCouplingFicksLaw
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FACET_COUPLING_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_FACET_COUPLING_FICKS_LAW_HH

#include <vector>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/flux/box/fickslaw.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Ficks's law for the box scheme in the context of coupled models
 *        where coupling occurs across the facets of the bulk domain elements
 *        with a lower-dimensional domain living on these facets.
 */
template<class TypeTag, ReferenceSystemFormulation referenceSystem = ReferenceSystemFormulation::massAveraged>
class BoxFacetCouplingFicksLaw
: public FicksLawImplementation<TypeTag, DiscretizationMethod::box, referenceSystem>
{
    using ParentType = FicksLawImplementation<TypeTag, DiscretizationMethod::box, referenceSystem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using BalanceEqOpts = GetPropType<TypeTag, Properties::BalanceEqOpts>;
    static constexpr int numComponents = ModelTraits::numFluidComponents();

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;

public:

    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static ComponentFluxVector flux(const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace& scvf,
                                    const int phaseIdx,
                                    const ElementFluxVarsCache& elemFluxVarCache)
    {
        // if this scvf is not on an interior boundary, use the standard law
        if (!scvf.interiorBoundary())
            return ParentType::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarCache);

        static const Scalar xi = getParamFromGroup<Scalar>(problem.paramGroup(), "FacetCoupling.Xi", 1.0);
        if ( !Dune::FloatCmp::eq(xi, 1.0, 1e-6) )
            DUNE_THROW(Dune::NotImplemented, "Xi != 1.0 cannot be used with the Box-Facet-Coupling scheme");

        // get some references for convenience
        const auto& fluxVarCache = elemFluxVarCache[scvf];
        const auto& shapeValues = fluxVarCache.shapeValues();
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        // interpolate density to scvf integration point
        Scalar rho = 0.0;
        for (const auto& scv : scvs(fvGeometry))
            rho += massOrMolarDensity(elemVolVars[scv], referenceSystem, phaseIdx)*shapeValues[scv.indexInElement()][0];

        // on interior Neumann boundaries, evaluate the flux using the facet effective diffusion coefficient
        ComponentFluxVector componentFlux(0.0);
        const auto bcTypes = problem.interiorBoundaryTypes(element, scvf);
        if (bcTypes.hasOnlyNeumann())
        {
            // compute tpfa flux from integration point to facet centerline
            const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);

            using std::sqrt;
            // If this is a surface grid, use the square root of the facet extrusion factor
            // as an approximate average distance from scvf ip to facet center
            using std::sqrt;
            const auto a = facetVolVars.extrusionFactor();
            auto preGradX = scvf.unitOuterNormal();
            preGradX *= dim == dimWorld ? 0.5*a : 0.5*sqrt(a);
            preGradX /= preGradX.two_norm2();

            for (int compIdx = 0; compIdx < numComponents; compIdx++)
            {
                if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

                // interpolate mole fraction to scvf integration point
                Scalar x = 0.0;
                for (const auto& scv : scvs(fvGeometry))
                    x += massOrMoleFraction(elemVolVars[scv], referenceSystem, phaseIdx, compIdx)*shapeValues[scv.indexInElement()][0];

                // compute the diffusive flux by means of a finite difference
                auto gradX = preGradX;
                gradX *= (massOrMoleFraction(facetVolVars, referenceSystem, phaseIdx, compIdx) - x);

                componentFlux[compIdx] = -1.0*rho*Extrusion::area(scvf)
                                             *insideVolVars.extrusionFactor()
                                             *vtmv(scvf.unitOuterNormal(),
                                                   facetVolVars.effectiveDiffusionCoefficient(phaseIdx, phaseIdx, compIdx),
                                                   gradX);

                if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx) && !FluidSystem::isTracerFluidSystem())
                    componentFlux[FluidSystem::getMainComponent(phaseIdx)] -= componentFlux[compIdx];
            }

            return componentFlux;
        }

        // on interior Dirichlet boundaries use the facet mass/mole fraction and evaluate flux
        else if (bcTypes.hasOnlyDirichlet())
        {
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
            {
                if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

                // create vector with nodal mole/mass fractions
                std::vector<Scalar> xFractions(element.subEntities(dim));
                for (const auto& scv : scvs(fvGeometry))
                    xFractions[scv.localDofIndex()] = massOrMoleFraction(elemVolVars[scv], referenceSystem, phaseIdx, compIdx);

                // substitute with facet mole/mass fractions for those scvs touching this facet
                for (const auto& scvfJ : scvfs(fvGeometry))
                    if (scvfJ.interiorBoundary() && scvfJ.facetIndexInElement() == scvf.facetIndexInElement())
                        xFractions[ fvGeometry.scv(scvfJ.insideScvIdx()).localDofIndex() ]
                                 = massOrMoleFraction(problem.couplingManager().getLowDimVolVars(element, scvfJ), referenceSystem, phaseIdx, compIdx);

                // evaluate gradX at integration point
                Dune::FieldVector<Scalar, dimWorld> gradX(0.0);
                for (const auto& scv : scvs(fvGeometry))
                    gradX.axpy(xFractions[scv.localDofIndex()], fluxVarCache.gradN(scv.indexInElement()));

                // apply matrix diffusion coefficient and return the flux
                componentFlux[compIdx] = -1.0*rho*Extrusion::area(scvf)
                                             *insideVolVars.extrusionFactor()
                                             *vtmv(scvf.unitOuterNormal(),
                                                   insideVolVars.effectiveDiffusionCoefficient(phaseIdx, phaseIdx, compIdx),
                                                   gradX);

                if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx) && !FluidSystem::isTracerFluidSystem())
                    componentFlux[FluidSystem::getMainComponent(phaseIdx)] -= componentFlux[compIdx];
            }

            return componentFlux;
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
        DUNE_THROW(Dune::NotImplemented, "transmissibilty computation for BoxFacetCouplingFicksLaw");
    }
};

} // end namespace Dumux

#endif
