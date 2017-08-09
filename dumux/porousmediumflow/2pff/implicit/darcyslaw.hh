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
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_FRACFLOW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_FRACFLOW_HH

#include <memory>

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/darcyslaw.hh>

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
class FractionalFlowDarcysLaw
{
    using Implementation = DarcysLawImplementation<TypeTag, DiscretizationMethods::CCTpfa>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    enum {
        viscousFluxIdx = 0,
        gravityFluxIdx = 1,
        capillaryFluxIdx = 2
    };

    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCTpfa;

    // state the type for the corresponding cache and its filler
    using Cache = typename Implementation::Cache;
    using CacheFiller = typename Implementation::CacheFiller;

    // the return type of the flux method
    // we return three separate fluxes: the viscous, gravity and capillary flux
    using ReturnType = std::array<Scalar, 3>;

    static ReturnType flux(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf,
                           int phaseIdx,
                           const ElementFluxVarsCache& elemFluxVarsCache)
    {
        ReturnType fluxes;

        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        // Get the total velocity from the input file
        static const GlobalPosition v_t = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, Problem, TotalVelocity);

        // The viscous flux before upwinding is the total velocity
        fluxes[viscousFluxIdx] = (v_t*scvf.unitOuterNormal())*scvf.area(); // TODO extrusion factor

        // Compute the capillary flux before upwinding
        const auto sInside = insideVolVars.saturation(wPhaseIdx);
        const auto sOutside = outsideVolVars.saturation(wPhaseIdx);
        fluxes[capillaryFluxIdx] = fluxVarsCache.tij()*(sInside - sOutside);


        // Compute the gravitational flux before upwinding
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
        {
            // do averaging for the density over all neighboring elements
            const auto rho = [&](int phaseIdx)
            {
                // boundaries
                if (scvf.boundary())
                    return insideVolVars.density(phaseIdx);

                // inner faces with two neighboring elements
                else
                    return (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;
            };

            // ask for the gravitational acceleration in the inside neighbor
            const auto xInside = insideScv.center();
            const auto gInside = problem.gravityAtPos(xInside);
            const auto xOutside = scvf.boundary() ? scvf.ipGlobal()
                                                  : fvGeometry.scv(scvf.outsideScvIdx()).center();
            const auto gOutside = problem.gravityAtPos(xOutside);

            fluxes[gravityFluxIdx] = fluxVarsCache.tij()*(xInside*gInside - xOutside*gOutside)*(rho(wPhaseIdx)- rho(nPhaseIdx));
        }
        else // no gravity
        {
            fluxes[gravityFluxIdx] = 0.0;
        }

        return fluxes;
    }
};

} // end namespace Dumux

#endif
