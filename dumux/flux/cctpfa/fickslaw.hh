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
 * \ingroup CCTpfaFlux
 * \brief Fick's law for cell-centered finite volume schemes with two-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FICKS_LAW_HH

#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/flux/fickiandiffusioncoefficients.hh>

#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation;

/*!
 * \ingroup CCTpfaFlux
 * \brief Fick's law for cell-centered finite volume schemes with two-point flux approximation
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation<TypeTag, DiscretizationMethod::cctpfa, referenceSystem>
{
    using Implementation = FicksLawImplementation<TypeTag, DiscretizationMethod::cctpfa, referenceSystem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BalanceEqOpts = GetPropType<TypeTag, Properties::BalanceEqOpts>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    static const int numPhases = ModelTraits::numFluidPhases();
    static const int numComponents = ModelTraits::numFluidComponents();

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
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;
    //! Return the reference system
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    //! state the type for the corresponding cache and its filler
    using Cache = TpfaFicksLawCache;

    using DiffusionCoefficientsContainer = FickianDiffusionCoefficients<Scalar, numPhases, numComponents>;

    /*!
     * \brief Returns the diffusive fluxes of all components within
     *        a fluid phase across the given sub-control volume face.
     *        The computed fluxes are given in mole/s or kg/s, depending
     *        on the template parameter ReferenceSystemFormulation.
     */
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
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

            // diffusion tensors are always solution dependent
            Scalar tij = elemFluxVarsCache[scvf].diffusionTij(phaseIdx, compIdx);

            // get inside/outside volume variables
            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

            // the inside and outside mass/mole fractions fractions
            const Scalar xInside = massOrMoleFraction(insideVolVars, referenceSystem, phaseIdx, compIdx);
            const Scalar massOrMoleFractionOutside = massOrMoleFraction(outsideVolVars, referenceSystem, phaseIdx, compIdx);
            const Scalar xOutside = scvf.numOutsideScvs() == 1 ? massOrMoleFractionOutside
                                  : branchingFacetX(problem, element, fvGeometry, elemVolVars,
                                                    elemFluxVarsCache, scvf, xInside, tij, phaseIdx, compIdx);

            const Scalar rhoInside = massOrMolarDensity(insideVolVars, referenceSystem, phaseIdx);
            const Scalar rhoOutside = massOrMolarDensity(outsideVolVars, referenceSystem, phaseIdx);

            const Scalar rho = scvf.numOutsideScvs() == 1 ? 0.5*(rhoInside + rhoOutside)
                                                        : branchingFacetDensity(elemVolVars, scvf, phaseIdx, rhoInside);

            componentFlux[compIdx] = rho*tij*(xInside - xOutside);
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx))
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


        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const auto getDiffCoeff = [&](const auto& vv)
        {
            if constexpr (FluidSystem::isTracerFluidSystem())
                return vv.effectiveDiffusionCoefficient(0, 0, compIdx);
            else
                return vv.effectiveDiffusionCoefficient(phaseIdx, FluidSystem::getMainComponent(phaseIdx), compIdx);
        };

        const auto insideD = getDiffCoeff(insideVolVars);

        const Scalar ti = computeTpfaTransmissibility(scvf, insideScv, insideD, insideVolVars.extrusionFactor());

        // for the boundary (dirichlet) or at branching points we only need ti
        Scalar tij;
        if (scvf.boundary() || scvf.numOutsideScvs() > 1)
            tij = Extrusion::area(scvf)*ti;

        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const auto outsideD = getDiffCoeff(outsideVolVars);

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
                tij = Extrusion::area(scvf)*(ti * tj)/(ti + tj);
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
            const Scalar massOrMoleFractionOutside = massOrMoleFraction(outsideVolVars, referenceSystem, phaseIdx, compIdx);
            const auto& flippedScvf = fvGeometry.flipScvf(scvf.index(), i);

            const Scalar outsideTi = elemFluxVarsCache[flippedScvf].diffusionTij(phaseIdx, compIdx);
            sumTi += outsideTi;
            sumXTi += outsideTi*massOrMoleFractionOutside;
        }

        return sumTi > 0 ? sumXTi/sumTi : 0;
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
            const Scalar rhoOutside = massOrMolarDensity(outsideVolVars, referenceSystem, phaseIdx);
            rho += rhoOutside;
        }
        return rho/(scvf.numOutsideScvs()+1);
    }
};

} // end namespace Dumux

#endif
