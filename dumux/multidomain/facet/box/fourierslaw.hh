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
 * \copydoc Dumux::BoxFacetCouplingFouriersLaw
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FACET_COUPLING_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_FACET_COUPLING_FOURIERS_LAW_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/box/fourierslaw.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief Fourier's law for the box scheme scheme in the context of coupled models
 *        where coupling occurs across the facets of the bulk domain elements
 *        with a lower-dimensional domain living on these facets.
 *
 * \tparam TypeTag the problem type tag
 */
template<class TypeTag>
class BoxFacetCouplingFouriersLaw
{
    using DefaultBoxFouriersLaw = FouriersLawImplementation<TypeTag, DiscretizationMethod::box>;

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

public:

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
            return DefaultBoxFouriersLaw::flux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarCache);

        // obtain the effective diffusivity model used
        using EffThermCondModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);

        // get some references for convenience
        const auto& fluxVarCache = elemFluxVarCache[scvf];
        const auto& shapeValues = fluxVarCache.shapeValues();
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        // evaluate user-defined interior boundary types
        const auto bcTypes = problem.interiorBoundaryTypes(element, scvf);

        // on interior Neumann boundaries, evaluate the flux using the facet permeability
        if (bcTypes.hasOnlyNeumann())
        {
            // interpolate temperature to scvf integration point
            Scalar T = 0.0;
            for (const auto& scv : scvs(fvGeometry))
                T += elemVolVars[scv].temperature()*shapeValues[scv.indexInElement()][0];

            // compute tpfa flux from integration point to facet centerline
            const auto& facetData = problem.couplingManager().getLowDimCouplingData(element, scvf);

            using std::sqrt;
            // If this is a surface grid, use the square root of the facet extrusion factor
            // as an approximate average distance from scvf ip to facet center
            const auto a = facetData.volVars().extrusionFactor();
            auto gradT = scvf.unitOuterNormal();
            gradT *= dim == dimWorld ? 0.5*a : 0.5*sqrt(a);
            gradT /= gradT.two_norm2();
            gradT *= (facetData.volVars().temperature() - T);

            // TODO How can we tell this class which effective thermal conductivity law is used
            //      in the facet domain? Here, we simply assume it uses the same law
            //      as this domain. But, we cannot make it an additional template parameter
            //      because that leads to a compiler error for models that do not specify
            //      an effective diffusivity law, e.g. models that do not consider diffusion.
            const auto facetLambda = EffThermCondModel::effectiveThermalConductivity(facetData.volVars(),
                                                                                     facetData.problem().spatialParams(),
                                                                                     facetData.element(),
                                                                                     facetData.fvGeometry(),
                                                                                     facetData.scv());

            // apply facet fourier coefficient
            return -1.0*scvf.area()*insideVolVars.extrusionFactor()*vtmv(scvf.unitOuterNormal(), facetLambda, gradT);
        }

        // on interior Dirichlet boundaries use the facet mole fraction and evaluate flux
        else if (bcTypes.hasOnlyDirichlet())
        {
                // create vector with nodal temperatures
                std::vector<Scalar> T(element.subEntities(dim));
                for (const auto& scv : scvs(fvGeometry))
                    T[scv.localDofIndex()] = elemVolVars[scv].temperature();

                // substitute with facet pressures for those scvs touching this facet
                for (const auto& scvfJ : scvfs(fvGeometry))
                    if (scvfJ.interiorBoundary() && scvfJ.facetIndexInElement() == scvf.facetIndexInElement())
                        T[ fvGeometry.scv(scvfJ.insideScvIdx()).localDofIndex() ]
                                 = problem.couplingManager().getLowDimVolVars(element, scvfJ).temperature();

                // evaluate gradT at integration point
                Dune::FieldVector<Scalar, dimWorld> gradT(0.0);
                for (const auto& scv : scvs(fvGeometry))
                    gradT.axpy(T[scv.localDofIndex()], fluxVarCache.gradN(scv.indexInElement()));

                auto insideLambda = EffThermCondModel::effectiveThermalConductivity(insideVolVars,
                                                                                    problem.spatialParams(),
                                                                                    element,
                                                                                    fvGeometry,
                                                                                    fvGeometry.scv(scvf.insideScvIdx()));

                // apply matrix fourier coefficient and return the flux
                return -1.0*scvf.area()*insideVolVars.extrusionFactor()*vtmv(scvf.unitOuterNormal(), insideLambda, gradT);
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
                                                           const FluxVarCache& fluxVarCache)
    {
        DUNE_THROW(Dune::NotImplemented, "transmissibilty computation for BoxFacetCouplingFouriersLaw");
    }
};

} // end namespace Dumux

#endif
