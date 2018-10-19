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
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \copydoc Dumux::BoxFacetCouplingFicksLaw
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FACET_COUPLING_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_FACET_COUPLING_FICKS_LAW_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/box/fickslaw.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief Fick's law for the box scheme scheme in the context of coupled models
 *        where coupling occurs across the facets of the bulk domain elements
 *        with a lower-dimensional domain living on these facets.
 *
 * \tparam TypeTag the problem type tag
 */
template<class TypeTag>
class BoxFacetCouplingFicksLaw
{
    using DefaultBoxFicksLaw = FicksLawImplementation<TypeTag, DiscretizationMethod::box>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BalanceEqOpts = typename GET_PROP_TYPE(TypeTag, BalanceEqOpts);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using CoordScalar = typename GridView::ctype;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    static constexpr int numPhases = ModelTraits::numPhases();
    static constexpr int numComponents = ModelTraits::numComponents();
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
            return DefaultBoxFicksLaw::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarCache);

        // obtain the effective diffusivity model used
        using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);

        // get some references for convenience
        const auto& fluxVarCache = elemFluxVarCache[scvf];
        const auto& shapeValues = fluxVarCache.shapeValues();
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        // evaluate user-defined interior boundary types
        const auto bcTypes = problem.interiorBoundaryTypes(element, scvf);

        ComponentFluxVector componentFlux(0.0);
        // on interior Neumann boundaries, evaluate the flux using the facet permeability
        if (bcTypes.hasOnlyNeumann())
        {
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
            {
                if(compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

                // interpolate mole fraction/density to scvf integration point
                Scalar x = 0.0;
                Scalar rho = 0.0;
                for (const auto& scv : scvs(fvGeometry))
                {
                    const auto& volVars = elemVolVars[scv];
                    x += volVars.moleFraction(phaseIdx, compIdx)*shapeValues[scv.indexInElement()][0];
                    rho += volVars.molarDensity(phaseIdx)*shapeValues[scv.indexInElement()][0];
                }

                // compute tpfa flux from integration point to facet centerline
                const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);

                using std::sqrt;
                // If this is a surface grid, use the square root of the facet extrusion factor
                // as an approximate average distance from scvf ip to facet center
                const auto a = facetVolVars.extrusionFactor();
                auto gradX = scvf.unitOuterNormal();
                gradX *= dim == dimWorld ? 0.5*a : 0.5*sqrt(a);
                gradX /= gradX.two_norm2();
                gradX *= (facetVolVars.moleFraction(phaseIdx, compIdx) - x);

                // TODO How can we tell this class which effective diffusivity law is used
                //      in the facet domain? Here, we simply assume it uses the same law
                //      as this domain. But, we cannot make it an additional template parameter
                //      because that leads to a compiler error for models that do not specify
                //      an effective diffusivity law, e.g. models that do not consider diffusion.
                auto facetD = EffDiffModel::effectiveDiffusivity(facetVolVars.porosity(),
                                                                 facetVolVars.saturation(phaseIdx),
                                                                 facetVolVars.diffusionCoefficient(phaseIdx, compIdx));

                // apply facet diffusion coefficient
                componentFlux[compIdx] = -1.0*rho*scvf.area()
                                             *insideVolVars.extrusionFactor()
                                             *vtmv(scvf.unitOuterNormal(), facetD, gradX);
                if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx) && !FluidSystem::isTracerFluidSystem())
                    componentFlux[phaseIdx] -= componentFlux[compIdx];
            }
        }

        // on interior Dirichlet boundaries use the facet mole fraction and evaluate flux
        else if (bcTypes.hasOnlyDirichlet())
        {
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
            {
                if(compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

                // create vector with nodal mole fractions
                std::vector<Scalar> x(element.subEntities(dim));
                for (const auto& scv : scvs(fvGeometry))
                    x[scv.localDofIndex()] = elemVolVars[scv].moleFraction(phaseIdx, compIdx);

                // substitute with facet pressures for those scvs touching this facet
                for (const auto& scvfJ : scvfs(fvGeometry))
                    if (scvfJ.interiorBoundary() && scvfJ.facetIndexInElement() == scvf.facetIndexInElement())
                        x[ fvGeometry.scv(scvfJ.insideScvIdx()).localDofIndex() ]
                                 = problem.couplingManager().getLowDimVolVars(element, scvfJ).moleFraction(phaseIdx, compIdx);

                // evaluate gradX at integration point
                Scalar rho(0.0);
                Dune::FieldVector<Scalar, dimWorld> gradX(0.0);
                for (const auto& scv : scvs(fvGeometry))
                {
                    rho += elemVolVars[scv].molarDensity(phaseIdx)*shapeValues[scv.indexInElement()][0];
                    gradX.axpy(x[scv.localDofIndex()], fluxVarCache.gradN(scv.indexInElement()));
                }

                auto insideD = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(),
                                                                  insideVolVars.saturation(phaseIdx),
                                                                  insideVolVars.diffusionCoefficient(phaseIdx, compIdx));

                // apply matrix diffusion coefficient and return the flux
                componentFlux[compIdx] = -1.0*scvf.area()*rho
                                             *insideVolVars.extrusionFactor()
                                             *vtmv(scvf.unitOuterNormal(), insideD, gradX);
                if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx) && !FluidSystem::isTracerFluidSystem())
                    componentFlux[phaseIdx] -= componentFlux[compIdx];
            }
        }

        // mixed boundary types are not supported
        else
            DUNE_THROW(Dune::NotImplemented, "Mixed boundary types are not supported");

        // return computed fluxes
        return componentFlux;
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
        DUNE_THROW(Dune::NotImplemented, "transmissibilty computation for BoxFacetCouplingFicksLaw");
    }
};

} // end namespace Dumux

#endif
