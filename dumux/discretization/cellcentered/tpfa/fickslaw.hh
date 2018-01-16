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
 * \ingroup CCTpfaDiscretization
 * \brief Fick's law for cell-centered finite volume schemes with two-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FICKS_LAW_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

namespace Dumux
{
// forward declaration
template<class TypeTag, DiscretizationMethods discMethod>
class FicksLawImplementation;

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Fick's law for cell-centered finite volume schemes with two-point flux approximation
 */
template <class TypeTag>
class FicksLawImplementation<TypeTag, DiscretizationMethods::CCTpfa>
{
    using Implementation = FicksLawImplementation<TypeTag, DiscretizationMethods::CCTpfa>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BalanceEqOpts = typename GET_PROP_TYPE(TypeTag, BalanceEqOpts);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static const int numComponents = GET_PROP_VALUE(TypeTag,NumComponents);

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;

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

    //! Class that caches the transmissibility
    class TpfaFicksLawCache
    {
    public:
        using Filler = TpfaFicksLawCacheFiller;

        void updateDiffusion(const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace &scvf,
                             const unsigned int phaseIdx,
                             const unsigned int compIdx)
        {
            tij_[phaseIdx][compIdx] = Implementation::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, compIdx);
        }

        const Scalar& diffusionTij(unsigned int phaseIdx, unsigned int compIdx) const
        { return tij_[phaseIdx][compIdx]; }

    private:
        std::array< std::array<Scalar, numComponents>, numPhases> tij_;
    };

public:
    //! state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCTpfa;

    //! state the type for the corresponding cache and its filler
    using Cache = TpfaFicksLawCache;

    //! return diffusive fluxes for all components in a phase
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
            if(compIdx == FluidSystem::getMainComponent(phaseIdx))
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
            if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx) && !FluidSystem::isTracerFluidSystem())
                componentFlux[FluidSystem::getMainComponent(phaseIdx)] -= componentFlux[compIdx];
        }

        return componentFlux;
    }

    //! compute diffusive transmissibilities
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

        using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
        const auto insideD = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(),
                                                                insideVolVars.saturation(phaseIdx),
                                                                insideVolVars.diffusionCoefficient(phaseIdx, compIdx));
        const Scalar ti = computeTpfaTransmissibility(scvf, insideScv, insideD, insideVolVars.extrusionFactor());

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

            const auto outsideD = EffDiffModel::effectiveDiffusivity(outsideVolVars.porosity(),
                                                                     outsideVolVars.saturation(phaseIdx),
                                                                     outsideVolVars.diffusionCoefficient(phaseIdx, compIdx));

            Scalar tj;
            if (dim == dimWorld)
                // assume the normal vector from outside is anti parallel so we save flipping a vector
                tj = -1.0*computeTpfaTransmissibility(scvf, outsideScv, outsideD, outsideVolVars.extrusionFactor());
            else
                tj = computeTpfaTransmissibility(fvGeometry.flipScvf(scvf.index()),
                                                 outsideScv,
                                                 outsideD,
                                                 outsideVolVars.extrusionFactor());

            // check if we are dividing by zero!
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = scvf.area()*(ti * tj)/(ti + tj);
        }

        return tij;
    }

private:
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
};

} // end namespace Dumux

#endif
