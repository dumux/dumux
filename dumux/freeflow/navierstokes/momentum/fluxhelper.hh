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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMomentumBoundaryFluxHelper
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_BOUNDARY_FLUXHELPER_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_BOUNDARY_FLUXHELPER_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/momentum/velocitygradients.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Struct containing flux helper functions to be used in the momentum problem's Neumann function.
 */
struct NavierStokesMomentumBoundaryFluxHelper
{
    /*!
     * \brief Returns the momentum flux a fixed-pressure boundary.
     * \param problem The problem
     * \param fvGeometry The finite-volume geometry
     * \param scvf The sub control volume face
     * \param elemVolVars The volume variables for the element
     * \param elemFluxVarsCache The flux variables cache for the element
     * \param pressure The pressure given at the boundary
     * \param zeroNormalVelocityGradient If set to false, this yields a "natural" outflow condition where the
     *                                   flow tends to diverge at the boundary. If set to true (default), this
     *                                   yields a "do-nothing" outflow condition where \f$\mathbf{v} \cdot \mathbf{n} = 0\f$
     *                                   is satisfied and the general flow field is preserved in many cases.
     *                                   See, e.g., https://www.math.uni-magdeburg.de/~richter/WS17/numns/fem.pdf for details.
     */
    template<class Problem, class FVElementGeometry, class ElementVolumeVariables, class ElementFluxVariablesCache, class Scalar>
    static auto fixedPressureMomentumFlux(const Problem& problem,
                                          const FVElementGeometry& fvGeometry,
                                          const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                          const ElementVolumeVariables& elemVolVars,
                                          const ElementFluxVariablesCache& elemFluxVarsCache,
                                          const Scalar pressure,
                                          const bool zeroNormalVelocityGradient = true)
    {
        // TODO density upwinding?
        static_assert(FVElementGeometry::GridGeometry::discMethod == DiscretizationMethod::fcstaggered);
        using NumEqVector = typename Problem::Traits::NumEqVector;
        NumEqVector flux(0.0);
        const auto& element = fvGeometry.element();

        if (scvf.isFrontal())
        {
            // pressure contribution
            flux[scvf.normalAxis()] = (pressure - problem.referencePressure(element, fvGeometry, scvf)) * scvf.directionSign();

            if (problem.enableInertiaTerms())
            {
                const auto v = elemVolVars[scvf.insideScvIdx()].velocity();
                flux[scvf.normalAxis()] += v*v * problem.density(element, fvGeometry, scvf) * scvf.directionSign();
            }
        }

        // if no zero velocity gradient is desired we are done here, ....
        if (!zeroNormalVelocityGradient)
            return flux;

        // ..., otherwise, make sure the flow does not diverge by accounting for the off-diagonal entries of the stress tensor
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

        if (scvf.isLateral())
        {
            // viscous terms
            const Scalar mu = problem.effectiveViscosity(element, fvGeometry, scvf);

            // lateral face normal to boundary (integration point touches boundary)
            if (scv.boundary())
                flux[scv.dofAxis()] -= mu * StaggeredVelocityGradients::velocityGradIJ(fvGeometry, scvf, elemVolVars)
                                               * scvf.directionSign();
            // lateral face coinciding with boundary
            else if (scvf.boundary())
                flux[scv.dofAxis()] -= mu * StaggeredVelocityGradients::velocityGradJI(fvGeometry, scvf, elemVolVars)
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
                        if (const auto bcTypes = problem.boundaryTypes(element, scvf); bcTypes.isDirichlet(scvf.normalAxis()))
                            return problem.dirichlet(element, scvf)[scvf.normalAxis()];
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
                    const auto outerVelocity = [&]
                    {
                        if (scvf.boundary())
                        {
                            assert(problem.boundaryTypes(element, scvf).isDirichlet(scv.dofAxis()));
                            return problem.dirichlet(element, scvf)[scv.dofAxis()];
                        }
                        else
                            return elemVolVars[scvf.outsideScvIdx()].velocity();
                    }();

                    const auto rho = problem.getInsideAndOutsideDensity(element, fvGeometry, scvf);
                    const bool selfIsUpstream = scvf.directionSign() == sign(transportingVelocity);

                    const auto insideMomentum = innerVelocity * rho.first;
                    const auto outsideMomentum = outerVelocity * rho.second;

                    static const auto upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");

                    const auto transportedMomentum =  selfIsUpstream ? (upwindWeight * insideMomentum + (1.0 - upwindWeight) * outsideMomentum)
                                                                     : (upwindWeight * outsideMomentum + (1.0 - upwindWeight) * insideMomentum);

                    flux[scv.dofAxis()] += transportingVelocity * transportedMomentum * scvf.directionSign();
                }

                // lateral face coinciding with boundary
                else if (scvf.boundary())
                {
                    assert(false);
                    // this should not happen? TODO revise
                    // const auto insideDensity = problem.density(element, fvGeometry.scv(scvf.insideScvIdx()));
                    // const auto innerVelocity = elemVolVars[scvf.insideScvIdx()].velocity();
                    // flux[scv.dofAxis()] += innerVelocity * transportingVelocity * insideDensity * scvf.directionSign();
                }
            }
        }

        return flux;
    }

    /*!
     * \brief Returns the momentum flux for a boundary with a slip condition.
     * \param problem The problem
     * \param fvGeometry The finite-volume geometry
     * \param scvf The sub control volume face
     * \param elemVolVars The volume variables for the element
     * \param elemFluxVarsCache The flux variables cache for the element
     * \param tangentialVelocityGradient A user-specified velocity gradient in tangential direction of the boundary. Only used for certain situations
     *                                   where this gradient cannot be determined automatically (e.g., at lateral Neumann boundaries).
     */
    template<class Problem, class FVElementGeometry, class ElementVolumeVariables, class ElementFluxVariablesCache, class TangentialVelocityGradient = double>
    static auto slipVelocityMomentumFlux(const Problem& problem,
                                         const FVElementGeometry& fvGeometry,
                                         const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                         const ElementVolumeVariables& elemVolVars,
                                         const ElementFluxVariablesCache& elemFluxVarsCache,
                                         const TangentialVelocityGradient& tangentialVelocityGradient = TangentialVelocityGradient(0.0))
    {
        static_assert(FVElementGeometry::GridGeometry::discMethod == DiscretizationMethod::fcstaggered);
        using NumEqVector = typename Problem::Traits::NumEqVector;
        NumEqVector flux(0.0);

        // slip velocity momentum contribution only makes sense for lateral scvfs
        if (scvf.isFrontal())
            return flux;

        using Scalar = std::decay_t<decltype(flux[0])>;
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
        const auto& orthogonalScv = fvGeometry.scv(orthogonalScvf.insideScvIdx());

        if (scvf.boundary() && problem.onSlipBoundary(fvGeometry, fvGeometry.frontalScvfOnBoundary(orthogonalScv)))
        {
            /*
            *                     ---------#######:::::::::     x dof position   ***** porous boundary at bottom
            *       slip          |      ||      |       ::
            *       gradient      |      ||      |  v_i  ::     -- element
            *       ------->      |      || scv  x~~~~>  ::
            *       ------>       |      ||      |       ::     O position at which gradient is evaluated (integration point)
            *       ----->        |  v_j ^       |       :^
            *       ---->         -------|-######O::::::::| <----This velocity dof (outer v_j) does not exist if the scv itself lies on a
            *                     ***************~~~>******     non-Dirichlet boundary. In that case, use the given tangentialVelocityGradient.
            *                      frontal scvf  v_slip
            *                      on porous                    || and # staggered half-control-volume (own element)
            *                      boundary
            *                                                   :: staggered half-control-volume (neighbor element)
            *
            */
            const Scalar v_i = elemVolVars[scvf.insideScvIdx()].velocity();

            // viscous terms
            const Scalar mu = problem.effectiveViscosity(fvGeometry.element(), fvGeometry, scvf);
            const Scalar distance = (scv.dofPosition()- scvf.ipGlobal()).two_norm();
            const Scalar velocityGrad_ji = [&]
            {
                if (elemVolVars.hasVolVars(orthogonalScvf.outsideScvIdx()))
                    return StaggeredVelocityGradients::velocityGradJI(fvGeometry, scvf, elemVolVars);
                else
                    return tangentialVelocityGradient;
            }();

            const Scalar slipVelocity = problem.beaversJosephVelocity(fvGeometry, scvf, elemVolVars, velocityGrad_ji)[scv.dofAxis()]; // TODO rename to slipVelocity
            const Scalar velocityGrad_ij = (slipVelocity - v_i) / distance * scvf.directionSign();

            flux[scv.dofAxis()] -= mu * velocityGrad_ij*scvf.directionSign();
            flux[scv.dofAxis()] -= mu * velocityGrad_ji*scvf.directionSign();

            // advective terms
            if (problem.enableInertiaTerms())
            {
                // transporting velocity corresponds to v_j
                const auto transportingVelocity = [&]()
                {
                    const auto innerTransportingVelocity = elemVolVars[orthogonalScvf.insideScvIdx()].velocity();

                    if (!elemVolVars.hasVolVars(orthogonalScvf.outsideScvIdx()))
                        return innerTransportingVelocity; // fallback

                    const auto outerTransportingVelocity = elemVolVars[orthogonalScvf.outsideScvIdx()].velocity();

                    // if the orthogonal scvf lies on a boundary and if there are outside volvars, we assume that these come from a Dirichlet condition
                    if (orthogonalScvf.boundary())
                        return outerTransportingVelocity;

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
                }();

                const auto innerVelocity = elemVolVars[scvf.insideScvIdx()].velocity();
                const auto outerVelocity = slipVelocity;

                const auto rho = problem.getInsideAndOutsideDensity(fvGeometry.element(), fvGeometry, scvf);
                const bool selfIsUpstream = scvf.directionSign() == sign(transportingVelocity);

                const auto insideMomentum = innerVelocity * rho.first;
                const auto outsideMomentum = outerVelocity * rho.second;

                static const auto upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");

                const auto transportedMomentum =  selfIsUpstream ? (upwindWeight * insideMomentum + (1.0 - upwindWeight) * outsideMomentum)
                                                                    : (upwindWeight * outsideMomentum + (1.0 - upwindWeight) * insideMomentum);

                flux[scv.dofAxis()] += transportingVelocity * transportedMomentum * scvf.directionSign();
            }
        }
        else if (scv.boundary() && problem.onSlipBoundary(fvGeometry, fvGeometry.frontalScvfOnBoundary(scv)))
        {
            /*                                    *
            *                     ---------#######*           x dof position   ***** porous boundary at right side
            *                     |      ||      |*
            *                     |      ||      |*  v_i      -- element
            *                     |      || scv  x~~~~>
            *                     |      ||      |*           O position at which gradient is evaluated (integration point)
            *                     |  v_j ^       |*^
            *                     -------|-######O*| v_slip   || and # staggered half-control-volume (own element)
            *                     |              |*
            *                     |  neighbor    |*
            *                     |   element     ~~~~> v_i outside (Does not exist if lower later lateral scvf lies on non-Dirichlet
            *                     |              |*                  boundary. In that case, use the given tangentialVelocityGradient.)
            *                     |              |*
            *                     ----------------*           :: staggered half-control-volume (neighbor element)
            *                            ^
            *                            |   ^
            *        slip gradient       |   |   ^
            *                            |   |   |
            *
            */
            const Scalar v_j = elemVolVars[orthogonalScvf.insideScvIdx()].velocity();

            // viscous terms
            const Scalar mu = problem.effectiveViscosity(fvGeometry.element(), fvGeometry, scvf);
            const Scalar distance = (fvGeometry.scv(orthogonalScvf.insideScvIdx()).dofPosition()- scvf.ipGlobal()).two_norm();

            const Scalar velocityGrad_ij = [&]
            {
                if (elemVolVars.hasVolVars(scvf.outsideScvIdx()))
                    return StaggeredVelocityGradients::velocityGradIJ(fvGeometry, scvf, elemVolVars);
                else
                    return tangentialVelocityGradient;
            }();

            const Scalar slipVelocity = problem.beaversJosephVelocity(fvGeometry, orthogonalScvf, elemVolVars, velocityGrad_ij)[scvf.normalAxis()]; // TODO rename to slipVelocity
            const Scalar velocityGrad_ji = (slipVelocity - v_j) / distance * orthogonalScvf.directionSign();

            flux[scv.dofAxis()] -= mu * velocityGrad_ij*scvf.directionSign();
            flux[scv.dofAxis()] -= mu * velocityGrad_ji*scvf.directionSign();

            // advective terms
            if (problem.enableInertiaTerms())
            {
                // transporting velocity corresponds to v_j
                const auto transportingVelocity = slipVelocity;

                // if the scvf lies on a boundary and if there are outside volvars, we assume that these come from a Dirichlet condition
                if (scvf.boundary() && elemVolVars.hasVolVars(scvf.outsideScvIdx()))
                {
                    flux[scv.dofAxis()] += problem.density(fvGeometry.element(), scv)
                                        * elemVolVars[scvf.outsideScvIdx()].velocity() * transportingVelocity * scvf.directionSign(); // TODO revise density
                    return flux;
                }

                const auto innerVelocity = elemVolVars[scvf.insideScvIdx()].velocity();
                const auto outerVelocity = [&]
                {
                    if (!elemVolVars.hasVolVars(scvf.outsideScvIdx()))
                        return innerVelocity; // fallback
                    else
                        return elemVolVars[scvf.outsideScvIdx()].velocity();
                }();

                const auto rho = problem.getInsideAndOutsideDensity(fvGeometry.element(), fvGeometry, scvf);
                const bool selfIsUpstream = scvf.directionSign() == sign(transportingVelocity);

                const auto insideMomentum = innerVelocity * rho.first;
                const auto outsideMomentum = outerVelocity * rho.second;

                static const auto upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");

                const auto transportedMomentum =  selfIsUpstream ? (upwindWeight * insideMomentum + (1.0 - upwindWeight) * outsideMomentum)
                                                                    : (upwindWeight * outsideMomentum + (1.0 - upwindWeight) * insideMomentum);

                flux[scv.dofAxis()] += transportingVelocity * transportedMomentum * scvf.directionSign();
            }
        }

        return flux;
    }
};

} // end namespace Dumux

#endif