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

#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_BOUNDARY_FLUXHELPER_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_BOUNDARY_FLUXHELPER_HH


#include <mutex>
#include <type_traits>
#include <dune/common/float_cmp.hh>
#include <dune/common/std/type_traits.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/momentum/velocitygradients.hh>

namespace Dumux {


struct NavierStokesMomentumBoundaryFluxHelper
{
    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class ElementFluxVariablesCache, class Scalar>
    static auto fixedPressureMomentumFlux(const Problem& problem,
                                          const Element& element,
                                          const FVElementGeometry& fvGeometry,
                                          const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                          const ElementVolumeVariables& elemVolVars,
                                          const ElementFluxVariablesCache& elemFluxVarsCache,
                                          const Scalar pressure,
                                          const bool zeroNormalVelocityGradient = true)
    {
        // TODO density upwinding?
        static_assert(FVElementGeometry::GridGeometry::discMethod == DiscretizationMethod::fcstaggered);
        using NumEqVector = decltype(problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf));
        NumEqVector flux(0.0);

        if (scvf.isFrontal())
        {
            // pressure contribution
            flux[scvf.directionIndex()] = (pressure - problem.referencePressure(element, fvGeometry, scvf)) * scvf.directionSign();

            if (problem.enableInertiaTerms())
            {
                const auto v = elemVolVars[scvf.insideScvIdx()].velocity();
                flux[scvf.directionIndex()] += v*v * problem.density(element, fvGeometry, scvf) * scvf.directionSign();
            }
        }

        // if no zero velocity gradient is desired we are done here, ....
        if (!zeroNormalVelocityGradient)
            return flux;

        // ..., otherwise, make sure the flow does not diverge by accounting for the off-diagonal entries of the stress tensor
        // TODO normal viscous term?

        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

        if (scvf.isLateral())
        {

            // lateral face normal to boundary (integration point touches boundary)
            if (scv.boundary() && scvf.boundary())
            {
                // const auto& frontalScvfOnBoundary = fvGeometry.frontalScvfOnBoundary(scv);
                const auto bcTypes = problem.boundaryTypes(element, scvf);
                if (bcTypes.isNeumann(scv.directionIndex()))
                {
                    static std::mutex recursionPreventionMutex;
                    if (!recursionPreventionMutex.try_lock())
                        DUNE_THROW(Dune::InvalidStateException, "fixedPressureMomentumFlux() was called recursively. "\
                                   << "To prevent a stack-overflow error, the simulation is aborted. "\
                                   << "Double check your neumann() function to avoid ambiguities in your domain corners.");

                    const auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
                    flux[scv.directionIndex()] = neumannFluxes[scv.directionIndex()];
                    recursionPreventionMutex.unlock();

                    return flux;
                }
            }

            // viscous terms
            const Scalar mu = problem.effectiveViscosity(element, fvGeometry, scvf);

            // lateral face normal to boundary (integration point touches boundary)
            if (scv.boundary())
                flux[scvf.directionIndex()] -= mu * StaggeredVelocityGradients::velocityGradIJ(fvGeometry, scvf, elemVolVars)
                                               * scvf.directionSign();
            // lateral face coinciding with boundary
            else if (scvf.boundary())
                flux[scv.directionIndex()] -= mu * StaggeredVelocityGradients::velocityGradJI(fvGeometry, scvf, elemVolVars)
                                              * scvf.directionSign();

            // advective terms
            if (problem.enableInertiaTerms())
            {
                const auto transportingVelocity = [&]()
                {
                    const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
                    const auto innerTransportingVelocity = elemVolVars[orthogonalScvf.insideScvIdx()].velocity();

                    if (scvf.boundary())
                    {
                        if (const auto bcTypes = problem.boundaryTypes(element, scvf); bcTypes.isDirichlet(scvf.directionIndex()))
                            return problem.dirichlet(element, scvf)[scvf.directionIndex()];
                        else
                            return
                                innerTransportingVelocity; // fallback
                    }
                    else
                    {
                        static const bool useOldScheme = getParam<bool>("FreeFlow.UseOldTransportingVelocity", true); // TODO how to deprecate?
                        if (useOldScheme)
                            return innerTransportingVelocity;
                        else
                        {
                            // average the transporting velocity by weighting with the scv volumes
                            const auto insideVolume = fvGeometry.scv(orthogonalScvf.insideScvIdx()).volume();
                            const auto outsideVolume = fvGeometry.scv(orthogonalScvf.outsideScvIdx()).volume();
                            const auto outerTransportingVelocity = elemVolVars[orthogonalScvf.outsideScvIdx()].velocity();
                            return (insideVolume*innerTransportingVelocity + outsideVolume*outerTransportingVelocity) / (insideVolume + outsideVolume);
                        }
                    }
                }();

                // lateral face normal to boundary (integration point touches boundary)
                if (scv.boundary())
                {
                    const auto innerVelocity = elemVolVars[scvf.insideScvIdx()].velocity();
                    const auto outerVelocity = elemVolVars[scvf.outsideScvIdx()].velocity();
                    const auto rho = problem.getInsideAndOutsideDensity(element, fvGeometry, scvf);

                    const bool selfIsUpstream = scvf.directionSign() == sign(transportingVelocity);

                    const auto insideMomentum = innerVelocity * rho.first;
                    const auto outsideMomentum = outerVelocity * rho.second;

                    static const auto upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");

                    const auto transportedMomentum =  selfIsUpstream ? (upwindWeight * insideMomentum + (1.0 - upwindWeight) * outsideMomentum)
                                                                     : (upwindWeight * outsideMomentum + (1.0 - upwindWeight) * insideMomentum);

                    flux[scvf.directionIndex()] += transportingVelocity * transportedMomentum * scvf.directionSign();
                }

                // lateral face coinciding with boundary
                else if (scvf.boundary())
                {
                    const auto insideDensity = problem.density(element, fvGeometry.scv(scvf.insideScvIdx()));
                    const auto innerVelocity = elemVolVars[scvf.insideScvIdx()].velocity();
                    flux[scv.directionIndex()] += innerVelocity * transportingVelocity * insideDensity * scvf.directionSign();
                }
            }
        }

        return flux;
    }
};

} // end namespace Dumux

#endif
