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
 * \brief This file contains the data which is required to calculate volume and mass fluxes of
 *        fluid phases over a face of a finite volume by means of the Darcy approximation in the
 *        presence of lower dimensional (coupled) elements living on the this domain's element facets.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FACET_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FACET_DARCYS_LAW_HH

#include <dumux/discretization/cellcentered/mpfa/darcyslaw.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>

namespace Dumux
{

/*!
 * \ingroup DarcysLaw
 * \brief Specialization of Darcy's Law for the CCMpfa method with lower dimensional
 *        elements living on the bulk elements' facets.
 */
template <class TypeTag>
class CCMpfaFacetCouplingDarcysLaw : public DarcysLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    // Always use the dynamic type for vectors (compatibility with the boundary)
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using CoefficientVector = typename BoundaryInteractionVolume::Traits::Vector;

public:
    static Scalar interpolateDensity(const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& elemVolVars,
                                     const SubControlVolumeFace& scvf,
                                     const FluxVariablesCache& fluxVarsCache,
                                     const unsigned int phaseIdx)
    {
        static const bool gravity = GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity);

        if (!gravity)
            return Scalar(0.0);
        else
        {
            // Return facet density on interior boundaries
            if (fluxVarsCache.isInteriorBoundary())
                return fluxVarsCache.interiorBoundaryDataSelf().facetVolVars(fvGeometry).density(phaseIdx);

            // use arithmetic mean of the densities around the scvf
            if (!scvf.boundary())
            {
                Scalar rho = elemVolVars[scvf.insideScvIdx()].density(phaseIdx);
                for (auto outsideIdx : scvf.outsideScvIndices())
                    rho += elemVolVars[outsideIdx].density(phaseIdx);
                return rho/(scvf.outsideScvIndices().size()+1);
            }
            else
                return elemVolVars[scvf.outsideScvIdx()].density(phaseIdx);
        }
    }

    static Scalar computeInteriorBoundaryContribution(const Problem& problem,
                                                      const FVElementGeometry& fvGeometry,
                                                      const ElementVolumeVariables& elemVolVars,
                                                      const FluxVariablesCache& fluxVarsCache,
                                                      unsigned int phaseIdx, Scalar rho)
    {
        static const bool gravity = GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity);

        // obtain the transmissibilites associated with all pressures
        const auto& tij = fluxVarsCache.advectionTij();

        // the interior dirichlet boundaries local indices start after
        // the cell and the domain Dirichlet boundary pressures
        const auto startIdx = fluxVarsCache.advectionVolVarsStencil().size();

        // The vector of interior neumann fluxes
        const auto& cij = fluxVarsCache.advectionCij();
        CoefficientVector facetCouplingFluxes(cij.size(), 0.0);

        // add interior Dirichlet boundary contributions
        Scalar flux = 0.0;
        for (auto&& data : fluxVarsCache.interiorBoundaryData())
        {
            // add neumann contributions
            if (data.faceType() == MpfaFaceTypes::interiorNeumann)
            {
                // get the scvf corresponding to actual interior neumann face
                const auto& curScvf = fvGeometry.scvf(data.scvfIndex());

                // get the volvars of the actual interior neumann face
                const auto facetVolVars = data.facetVolVars(fvGeometry, curScvf);

                // calculate "leakage factor"
                const auto n = curScvf.unitOuterNormal();
                const auto v = [&] ()
                                {
                                    auto res = n;
                                    res *= -0.5*facetVolVars.extrusionFactor();
                                    res /= res.two_norm2();
                                    return res;
                                } ();

                const auto facetH = [&] ()
                                    {
                                        if (!gravity)
                                            return facetVolVars.pressure(phaseIdx);
                                        else
                                        {
                                            const auto x = fvGeometry.scvf(data.scvfIndex()).ipGlobal();
                                            const auto g = problem.gravityAtPos(x);

                                            return facetVolVars.pressure(phaseIdx) - rho*(g*x);
                                        }
                                    } ();

                // add value to vector of interior neumann fluxes
                facetCouplingFluxes[data.localIndexInInteractionVolume()] += facetH*
                                                                             curScvf.area()*
                                                                             elemVolVars[curScvf.insideScvIdx()].extrusionFactor()*
                                                                             MpfaHelper::nT_M_v(n, facetVolVars.permeability(), v);
            }
            // Add additional Dirichlet fluxes for interior Dirichlet faces
            else if (data.faceType() == MpfaFaceTypes::interiorDirichlet)
            {
                Scalar h = data.facetVolVars(fvGeometry).pressure(phaseIdx);

                if (gravity)
                {
                    const auto x = fvGeometry.scvf(data.scvfIndex()).ipGlobal();
                    const auto g = problem.gravityAtPos(x);

                    h -= rho*(g*x);
                }

                flux += tij[startIdx + data.localIndexInInteractionVolume()]*h;
            }
        }

        return flux + cij*facetCouplingFluxes;
    }
};

} // end namespace

#endif
