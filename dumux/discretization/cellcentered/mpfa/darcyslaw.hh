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
 *        volume and mass fluxes of fluid phases over a face of a finite volume by means
 *        of the Darcy approximation. Specializations are provided for the different discretization methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_DARCYS_LAW_HH

#include <memory>

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>


namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(ProblemEnableGravity);
}

/*!
 * \ingroup DarcysLaw
 * \brief Specialization of Darcy's Law for the CCTpfa method.
 */
template <class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVarsCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

public:

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const unsigned int phaseIdx,
                       const FluxVarsCache& fluxVarsCache)
    {
        const auto& volVarsStencil = fluxVarsCache.advectionVolVarsStencil(phaseIdx);
        const auto& volVarsPositions = fluxVarsCache.advectionVolVarsPositions(phaseIdx);
        const auto& tij = fluxVarsCache.advectionTij(phaseIdx);

        // calculate Tij*pj
        Scalar flux(0.0);
        unsigned int localIdx = 0;
        for (const auto volVarIdx : volVarsStencil)
        {
            const auto& volVars = elemVolVars[volVarIdx];
            Scalar h = volVars.pressure(phaseIdx);

            // if gravity is enabled, add gravitational acceleration
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {
                // interface density as arithmetic mean of the two neighbors
                Scalar rho;
                if (!scvf.boundary())
                {
                    rho = elemVolVars[scvf.outsideScvIdx()].density(phaseIdx);
                    rho += elemVolVars[scvf.insideScvIdx()].density(phaseIdx);
                    rho /= 2.0;
                }
                else
                    rho = elemVolVars[scvf.outsideScvIdx()].density(phaseIdx);

                // gravitational acceleration in the center of the actual element
                const auto x = volVarsPositions[localIdx];
                const auto g = problem.gravityAtPos(x);

                h -= rho*(g*x);
            }

            flux += tij[localIdx++]*h;
        }

        return flux + fluxVarsCache.advectionNeumannFlux(phaseIdx);
    }

    static Stencil stencil(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvf)
    {
        Stencil stencil;
        if(problem.model().globalFvGeometry().scvfTouchesBoundary(scvf))
        {
            const auto& ivSeed = problem.model().globalFvGeometry().boundaryInteractionVolumeSeed(scvf);
            const auto& localScvSeeds = ivSeed.scvSeeds();

            stencil.reserve(localScvSeeds.size());
            for (const auto& localScvSeed : ivSeed.scvSeeds())
                stencil.push_back(localScvSeed.globalIndex());
        }
        else
        {
            const auto& ivSeed = problem.model().globalFvGeometry().interactionVolumeSeed(scvf);
            const auto& localScvSeeds = ivSeed.scvSeeds();

            stencil.reserve(localScvSeeds.size());
            for (const auto& localScvSeed : ivSeed.scvSeeds())
                stencil.push_back(localScvSeed.globalIndex());
        }

        return stencil;
    }
};

} // end namespace

#endif
