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
#ifndef DUMUX_SHALLOW_WATER_ADVECTIVE_FLUX_HH
#define DUMUX_SHALLOW_WATER_ADVECTIVE_FLUX_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/methods.hh>
//#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

namespace Dumux
{
template<class TypeTag>
class ShallowWaterAdvectiveFluxCacheFiller
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

public:
    template<class FluxVariablesCacheFiller>
    static void fill(FluxVariablesCache& scvfFluxVarsCache,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolumeFace& scvf,
                     const FluxVariablesCacheFiller& fluxVarsCacheFiller)
    {}
};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief The cache corresponding to tpfa numerical flux
 */
template<class TypeTag>
class ShallowWaterAdvectiveFluxCache
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

public:
    using Filler = ShallowWaterAdvectiveFluxCacheFiller<TypeTag>;

    void updateAdvection(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace &scvf)
    {}
};

template<class TypeTag>
class ShallowWaterAdvectiveFlux
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

  public:
    //! state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCTpfa;

    //! state the type for the corresponding cache
    using Cache = ShallowWaterAdvectiveFluxCache<TypeTag>;

    //! Compute the advective flux
    static NumEqVector flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        NumEqVector flux(0.0);

//         static const bool enableGravity = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Problem.EnableGravity");
//
//         const auto& fluxVarsCache = elemFluxVarsCache[scvf];
//
//         // Get the inside and outside volume variables
//         const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
//         const auto& insideVolVars = elemVolVars[insideScv];
//         const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
//
//         if (enableGravity)
//         {
//             // do averaging for the density over all neighboring elements
//             const auto rho = scvf.boundary() ? outsideVolVars.density(phaseIdx)
//                                              : (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;
//
//             // Obtain inside and outside pressures
//             const auto pInside = insideVolVars.pressure(phaseIdx);
//             const auto pOutside = outsideVolVars.pressure(phaseIdx);
//
//             const auto& tij = fluxVarsCache.advectionTij();
//             const auto& g = problem.gravityAtPos(scvf.ipGlobal());
//
//             //! compute alpha := n^T*K*g
//             const auto alpha_inside = vtmv(scvf.unitOuterNormal(), insideVolVars.permeability(), g)*insideVolVars.extrusionFactor();
//
//             Scalar flux = tij*(pInside - pOutside) + rho*scvf.area()*alpha_inside;
//
//             //! On interior faces we have to add K-weighted gravitational contributions
//             if (!scvf.boundary())
//             {
//                 const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
//                 const auto outsideK = outsideVolVars.permeability();
//                 const auto outsideTi = computeTpfaTransmissibility(scvf, outsideScv, outsideK, outsideVolVars.extrusionFactor());
//                 const auto alpha_outside = vtmv(scvf.unitOuterNormal(), outsideK, g)*outsideVolVars.extrusionFactor();
//
//                 flux += rho*tij/outsideTi*(alpha_inside - alpha_outside);
//             }
//
//             return flux;
//         }
//         else
//         {
//             // Obtain inside and outside pressures
//             const auto pInside = insideVolVars.pressure(phaseIdx);
//             const auto pOutside = outsideVolVars.pressure(phaseIdx);
//
//             // return flux
//             return fluxVarsCache.advectionTij()*(pInside - pOutside);
//         }
        return flux;
    }
};

} // end namespace Dumux

#endif
