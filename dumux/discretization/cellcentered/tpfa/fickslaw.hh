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
 * \brief This file contains the data which is required to calculate
 *        diffusive mass fluxes due to molecular diffusion with Fick's law.
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FICKS_LAW_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fluxvariablescaching.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(EffectiveDiffusivityModel);
}

/*!
 * \ingroup CCTpfaFicksLaw
 * \brief Specialization of Fick's Law for the CCTpfa method.
 */
template <class TypeTag>
class FicksLawImplementation<TypeTag, DiscretizationMethods::CCTpfa >
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static const int numComponents = GET_PROP_VALUE(TypeTag,NumComponents);

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;

    class TpfaFicksLawCache
    {
    public:
        void updateDiffusion(const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace &scvf,
                             const unsigned int phaseIdx,
                             const unsigned int compIdx)
        {
            tij_[phaseIdx][compIdx] = calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, compIdx);
        }

        const Scalar& diffusionTij(unsigned int phaseIdx, unsigned int compIdx) const
        { return tij_[phaseIdx][compIdx]; }

    private:
        std::array< std::array<Scalar, numComponents>, numPhases> tij_;
    };

    //! Class that fills the cache corresponding to tpfa Fick's Law
    class TpfaFicksLawCacheFiller
    {
    public:
        //! Function to fill a TpfaFicksLawCache of a given scvf
        //! This interface has to be met by any diffusion-related cache filler class
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         unsigned int phaseIdx, unsigned int compIdx,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace& scvf,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
            scvfFluxVarsCache.updateDiffusion(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, compIdx);
        }
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCTpfa;

    //! state the type for the corresponding cache and its filler
    using Cache = TpfaFicksLawCache;
    using CacheFiller = TpfaFicksLawCacheFiller;

    static ComponentFluxVector flux(const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace& scvf,
                                    const int phaseIdx,
                                    const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        ComponentFluxVector componentFlux(0.0);
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            if(FluidSystem::isMainComponent(compIdx, phaseIdx))
                continue;

            // diffusion tensors are always solution dependent
            Scalar tij = elemFluxVarsCache[scvf].diffusionTij(phaseIdx, compIdx);

            // get inside/outside volume variables
            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

            // the inside and outside mole fractions
            const auto xInside = insideVolVars.moleFraction(phaseIdx, compIdx);
            const auto xOutside = scvf.numOutsideScvs() == 1 ? outsideVolVars.moleFraction(phaseIdx, compIdx)
                                : branchingFacetX(problem, element, fvGeometry, elemVolVars,
                                                   elemFluxVarsCache, scvf, xInside, tij, phaseIdx, compIdx);

            const auto rhoInside = insideVolVars.molarDensity(phaseIdx);
            const auto rho = scvf.numOutsideScvs() == 1 ? 0.5*(rhoInside + outsideVolVars.molarDensity(phaseIdx))
                                                        : branchingFacetDensity(elemVolVars, scvf, phaseIdx, rhoInside);

            componentFlux[compIdx] = rho*tij*(xInside - xOutside);
            if (!FluidSystem::isTracerFluidSystem())
                componentFlux[phaseIdx] -= componentFlux[compIdx];
        }
        return componentFlux ;
    }

    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const SubControlVolumeFace& scvf,
                                            const int phaseIdx, const int compIdx)
    {
        Scalar tij;

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];

        auto insideD = insideVolVars.diffusionCoefficient(phaseIdx, compIdx);
        insideD = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(), insideVolVars.saturation(phaseIdx), insideD);
        const Scalar ti = calculateOmega_(scvf,
                                    insideD,
                                    insideScv,
                                    insideVolVars.extrusionFactor());

        // for the boundary (dirichlet) or at branching points we only need ti
        if (scvf.boundary() || scvf.numOutsideScvs() > 1)
        {
            tij = scvf.area()*ti;
        }
        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];

            auto outsideD = outsideVolVars.diffusionCoefficient(phaseIdx, compIdx);
            outsideD = EffDiffModel::effectiveDiffusivity(outsideVolVars.porosity(), outsideVolVars.saturation(phaseIdx), outsideD);

            Scalar tj;
            if (dim == dimWorld)
                // assume the normal vector from outside is anti parallel so we save flipping a vector
                tj = -1.0*calculateOmega_(scvf,
                                          outsideD,
                                          outsideScv,
                                          outsideVolVars.extrusionFactor());
            else
                tj = calculateOmega_(fvGeometry.flipScvf(scvf.index()),
                                     outsideD,
                                     outsideScv,
                                     outsideVolVars.extrusionFactor());

            // check if we are dividing by zero!
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = scvf.area()*(ti * tj)/(ti + tj);
        }

        return tij;
    }

    //! compute the mole/mass fraction at branching facets for network grids
    static Scalar branchingFacetX(const Problem& problem,
                                  const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& elemVolVars,
                                  const ElementFluxVariablesCache& elemFluxVarsCache,
                                  const SubControlVolumeFace& scvf,
                                  const Scalar insideX, const Scalar insideTi,
                                  const int phaseIdx, const int compIdx)
    {
        Scalar sumTi(insideTi);
        Scalar sumXTi(insideTi*insideX);

        for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
        {
            const auto outsideScvIdx = scvf.outsideScvIdx(i);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const auto& flippedScvf = fvGeometry.flipScvf(scvf.index(), i);

            auto outsideTi = elemFluxVarsCache[flippedScvf].diffusionTij(phaseIdx, compIdx);
            sumTi += outsideTi;
            sumXTi += outsideTi*outsideVolVars.moleFraction(phaseIdx, compIdx);
        }
        return sumXTi/sumTi;
    }

    //! compute the density at branching facets for network grids as arithmetic mean
    static Scalar branchingFacetDensity(const ElementVolumeVariables& elemVolVars,
                                        const SubControlVolumeFace& scvf,
                                        const int phaseIdx,
                                        const Scalar insideRho)
    {
        Scalar rho(insideRho);
        for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
        {
            const auto outsideScvIdx = scvf.outsideScvIdx(i);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            rho += outsideVolVars.molarDensity(phaseIdx);
        }
        return rho/(scvf.numOutsideScvs()+1);
    }

private:

    static Scalar calculateOmega_(const SubControlVolumeFace& scvf,
                                  const DimWorldMatrix &D,
                                  const SubControlVolume &scv,
                                  const Scalar extrusionFactor)
    {
        GlobalPosition Dnormal;
        D.mv(scvf.unitOuterNormal(), Dnormal);

        auto distanceVector = scvf.ipGlobal();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = Dnormal * distanceVector;
        omega *= extrusionFactor;

        return omega;
    }

    static Scalar calculateOmega_(const SubControlVolumeFace& scvf,
                                  const Scalar D,
                                  const SubControlVolume &scv,
                                  const Scalar extrusionFactor)
    {
        auto distanceVector = scvf.ipGlobal();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = D * (distanceVector * scvf.unitOuterNormal());
        omega *= extrusionFactor;

        return omega;
    }
};
} // end namespace

#endif
