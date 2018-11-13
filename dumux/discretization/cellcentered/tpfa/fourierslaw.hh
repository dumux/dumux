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
* \brief Fourier's law for cell-centered finite volume schemes with two-point flux approximation
*/
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FOURIERS_LAW_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

namespace Dumux
{
// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class FouriersLawImplementation;

/*!
* \ingroup CCTpfaDiscretization
* \brief Fourier's law for cell-centered finite volume schemes with two-point flux approximation
*/
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethod::cctpfa>
{
    using Implementation = FouriersLawImplementation<TypeTag, DiscretizationMethod::cctpfa>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache)::LocalView;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);

    //! Class that fills the cache corresponding to tpfa Fick's Law
    class TpfaFouriersLawCacheFiller
    {
    public:
        //! Function to fill a TpfaFicksLawCache of a given scvf
        //! This interface has to be met by any diffusion-related cache filler class
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace& scvf,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
            scvfFluxVarsCache.updateHeatConduction(problem, element, fvGeometry, elemVolVars, scvf);
        }
    };

    //! Class that caches the transmissibility
    class TpfaFouriersLawCache
    {
    public:
        using Filler = TpfaFouriersLawCacheFiller;

        void updateHeatConduction(const Problem& problem,
                                  const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& elemVolVars,
                                  const SubControlVolumeFace &scvf)
        {
            tij_ = Implementation::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
        }

        const Scalar& heatConductionTij() const
        { return tij_; }

    private:
        Scalar tij_;
    };

public:
    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;

    //! export the type for the corresponding cache
    using Cache = TpfaFouriersLawCache;

    //! Compute the heat condution flux assuming thermal equilibrium
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        // heat conductivities are always solution dependent (?)
        Scalar tij = elemFluxVarsCache[scvf].heatConductionTij();

        // get the inside/outside temperatures
        const auto tInside = elemVolVars[scvf.insideScvIdx()].temperature();
        const auto tOutside = scvf.numOutsideScvs() == 1 ? elemVolVars[scvf.outsideScvIdx()].temperature()
                              : branchingFacetTemperature_(problem, element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf, tInside, tij);

        return tij*(tInside - tOutside);
    }

    //! Compute transmissibilities
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const SubControlVolumeFace& scvf)
    {
        Scalar tij;

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];

        const auto insideLambda = ThermalConductivityModel::effectiveThermalConductivity(insideVolVars,
                                                                                         problem.spatialParams(),
                                                                                         element,
                                                                                         fvGeometry,
                                                                                         insideScv);
        const Scalar ti = computeTpfaTransmissibility(scvf, insideScv, insideLambda, insideVolVars.extrusionFactor());

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
            const auto outsideElement = fvGeometry.fvGridGeometry().element(outsideScvIdx);

            const auto outsideLambda = ThermalConductivityModel::effectiveThermalConductivity(outsideVolVars,
                                                                                              problem.spatialParams(),
                                                                                              outsideElement,
                                                                                              fvGeometry,
                                                                                              outsideScv);
            Scalar tj;
            if (dim == dimWorld)
                // assume the normal vector from outside is anti parallel so we save flipping a vector
                tj = -1.0*computeTpfaTransmissibility(scvf, outsideScv, outsideLambda, outsideVolVars.extrusionFactor());
            else
                tj = computeTpfaTransmissibility(fvGeometry.flipScvf(scvf.index()), outsideScv, outsideLambda, outsideVolVars.extrusionFactor());

            // check for division by zero!
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = scvf.area()*(ti * tj)/(ti + tj);
        }

        return tij;
    }

private:

    //! compute the temperature at branching facets for network grids
    static Scalar branchingFacetTemperature_(const Problem& problem,
                                             const Element& element,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars,
                                             const ElementFluxVarsCache& elemFluxVarsCache,
                                             const SubControlVolumeFace& scvf,
                                             Scalar insideTemperature,
                                             Scalar insideTi)
    {
        Scalar sumTi(insideTi);
        Scalar sumTempTi(insideTi*insideTemperature);

        for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
        {
            const auto outsideScvIdx = scvf.outsideScvIdx(i);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const auto& flippedScvf = fvGeometry.flipScvf(scvf.index(), i);
            const auto& outsideFluxVarsCache = elemFluxVarsCache[flippedScvf];

            auto outsideTi = outsideFluxVarsCache.heatConductionTij();
            sumTi += outsideTi;
            sumTempTi += outsideTi*outsideVolVars.temperature();
        }
        return sumTempTi/sumTi;
    }
};

} // end namespace Dumux

#endif
