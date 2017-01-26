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
 *        heat conduction fluxes with Fourier's law for cell-centered MPFA models.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FOURIERS_LAW_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

/*!
 * \ingroup FouriersLaw
 * \brief Specialization of Fourier's Law for the CCMpfa method.
 */
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);

    // Always use the dynamic type for vectors (compatibility with the boundary)
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using DynamicVector = typename BoundaryInteractionVolume::Vector;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    static constexpr bool facetCoupling = GET_PROP_VALUE(TypeTag, MpfaFacetCoupling);
    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

    static constexpr int energyEqIdx = GET_PROP_TYPE(TypeTag, Indices)::energyEqIdx;

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCMpfa;

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& volVarsStencil = fluxVarsCache.heatConductionVolVarsStencil();
        const auto& tij = fluxVarsCache.heatConductionTij();

        const bool isInteriorBoundary = enableInteriorBoundaries && fluxVarsCache.isInteriorBoundary();
        // For interior Neumann boundaries when using Tpfa for Neumann boundary conditions, we simply
        // return the user-specified flux
        if (isInteriorBoundary
            && useTpfaBoundary
            && fluxVarsCache.interiorBoundaryDataSelf().faceType() == MpfaFaceTypes::interiorNeumann)
            return scvf.area()*
                   elemVolVars[scvf.insideScvIdx()].extrusionFactor()*
                   problem.neumann(element,
                                   fvGeometry,
                                   elemVolVars,
                                   scvf)[energyEqIdx];

        // calculate Tij*tj
        Scalar flux(0.0);
        unsigned int localIdx = 0;
        for (const auto volVarIdx : volVarsStencil)
            flux += tij[localIdx++]*elemVolVars[volVarIdx].temperature();

        // if no interior boundaries are present, return heat conduction flux
        if (!enableInteriorBoundaries)
            return useTpfaBoundary ? flux : flux + fluxVarsCache.heatNeumannFlux();

        //////////////////////////////////////////////////////////////////
        // Handle interior boundaries
        //////////////////////////////////////////////////////////////////

        // get coefficients to transform the vector of interior neumann boundary conditions
        const auto& cij = fluxVarsCache.heatConductionCij();

        // The Vector of interior neumann fluxes
        DynamicVector facetCouplingFluxes(cij.size(), 0.0);
        for (auto&& data : fluxVarsCache.interiorBoundaryData())
        {
            // Add additional Dirichlet fluxes for interior Dirichlet faces
            if (data.faceType() == MpfaFaceTypes::interiorDirichlet)
            {
                // The transmissibilities of interior dirichlet boundaries are placed at the end
                // So we simply keep incrementing the local index
                flux += tij[localIdx + data.localIndexInInteractionVolume()]*data.facetVolVars(fvGeometry).temperature();
            }

            // add neumann fluxes for interior Neumann faces
            if (data.faceType() == MpfaFaceTypes::interiorNeumann)
            {
                if (facetCoupling)
                {
                    // get the scvf corresponding to actual interior boundary face
                    const auto& curScvf = fvGeometry.scvf(data.scvfIndex());

                    // obtain the complete data on the facet element
                    const auto completeFacetData = data.completeCoupledFacetData();

                    // calculate "lekage factor"
                    const auto n = curScvf.unitOuterNormal();
                    const auto v = [&] ()
                                    {
                                        auto res = n;
                                        res *= -0.5*completeFacetData.volVars().extrusionFactor();
                                        res -= curScvf.ipGlobal();
                                        res += curScvf.facetCorner();
                                        res /= res.two_norm2();
                                        return res;
                                    } ();

                    // get the thermal conductivity in the facet element
                    const auto facetLambda = ThermalConductivityModel::effectiveThermalConductivity(completeFacetData.volVars(),
                                                                                                    completeFacetData.spatialParams(),
                                                                                                    completeFacetData.element(),
                                                                                                    completeFacetData.fvGeometry(),
                                                                                                    completeFacetData.scv());

                    // add value to vector of interior neumann fluxes
                    facetCouplingFluxes[data.localIndexInInteractionVolume()] += MpfaHelper::nT_M_v(n,
                                                                                                      facetLambda,
                                                                                                      v);
                }
            }
        }

        // return overall resulting flux
        const Scalar interiorNeumannFlux = facetCoupling ? cij*facetCouplingFluxes : 0.0;
        return useTpfaBoundary ?
               flux + interiorNeumannFlux :
               flux + interiorNeumannFlux + fluxVarsCache.heatNeumannFlux();
    }
};

} // end namespace Dumux

#endif
