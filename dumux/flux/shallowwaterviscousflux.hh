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
 * \ingroup Flux
 * \copydoc Dumux::ShallowWaterViscousMomentumFlux
 */
#ifndef DUMUX_FLUX_SHALLOW_WATER_VISCOUS_FLUX_HH
#define DUMUX_FLUX_SHALLOW_WATER_VISCOUS_FLUX_HH

#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief Computes the shallow water viscous momentum flux due to (turbulent) viscosity
 *        by adding all surrounding shear stresses.
 *        For now implemented strictly for 2D depth-averaged models (i.e. 3 equations)
 */
template<class PrimaryVariables, class NumEqVector, typename std::enable_if_t<NumEqVector::size() == 3, int> = 0>
class ShallowWaterViscousFlux
{

public:

    using Cache = FluxVariablesCaching::EmptyDiffusionCache;
    using CacheFiller = FluxVariablesCaching::EmptyCacheFiller;
    /*!
     * \ingroup Flux
     * \brief Compute the viscous momentum flux contribution from the interface
     *        shear stress
     *
     *        The viscous momentum flux
     *        \f[
     *        \int \int_{V} \mathbf{\nabla} \cdot \nu_t h \mathbf{\nabla} \mathbf{u} dV
     *        \f]
     *        is re-written using Gauss' divergence theorem to:
     *        \f[
     *        \int_{S_f} \nu_t h \mathbf{\nabla} \mathbf{u} \cdot \mathbf{n_f} dS
     *        \f]
     *
     * \todo The implementation now is the simplest one involving
     *       only direct neighbours. This implementation is not complete/
     *       correct on non-orthogonal meshes. A more complete implementation
     *       with a more elaborate stencil that also takes into account
     *       the non-orthogonal contributions can be considered at a later stage.
     */
    template<class Problem, class FVElementGeometry, class ElementVolumeVariables>
    static NumEqVector flux(const Problem& problem,
                            const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {
        using Scalar = typename PrimaryVariables::value_type;

        using std::sqrt;
        using std::max;

        NumEqVector localFlux(0.0);

        // Get the inside and outside volume variables
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const auto& unitNormal = scvf.unitOuterNormal();

        // Inside and outside subcontrolvolumes
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());

        // Distance vector between the two cell centers
        const auto& cellCenterToCellCenter = outsideScv.center() - insideScv.center();

        // Check whether the mixing-length turbulence model is used
        static const auto useMixingLengthTurbulenceModel = getParam<bool>("ShallowWater.UseMixingLengthTurbulenceModel", false);

        // Compute the average shear velocity on the face
        const Scalar ustar = [&](){
            // If the mixing-length turbulence model is used:
            // compute ustar based on the bottom friction
            if (useMixingLengthTurbulenceModel)
            {
                // Get the bottom shear stress in the two adjacent cells
                // Note that element is not needed in spatialParams().frictionLaw (should be removed). For now we simply pass the same element twice
                const auto& bottomShearStressInside = problem.spatialParams().frictionLaw(element, insideScv).shearStress(insideVolVars);
                const auto& bottomShearStressOutside = problem.spatialParams().frictionLaw(element, outsideScv).shearStress(outsideVolVars);
                const auto bottomShearStressInsideMag = bottomShearStressInside.two_norm();
                const auto bottomShearStressOutsideMag = bottomShearStressOutside.two_norm();

                // Use a harmonic average of the viscosity and the depth at the interface.
                const auto averageBottomShearStress = 2.0*(bottomShearStressInsideMag*bottomShearStressOutsideMag)/max(1.0e-8,bottomShearStressInsideMag + bottomShearStressOutsideMag);

                // Shear velocity possibly needed in mixing-length turbulence model in the computation of the diffusion flux
                return sqrt(averageBottomShearStress);
            }
            else
                return 0.0;
        }();

        // The left (inside) and right (outside) states
        const auto waterDepthLeft  = insideVolVars.waterDepth();
        const auto velocityXLeft   = insideVolVars.velocity(0);
        const auto velocityYLeft   = insideVolVars.velocity(1);
        const auto waterDepthRight = outsideVolVars.waterDepth();
        const auto velocityXRight  = outsideVolVars.velocity(0);
        const auto velocityYRight  = outsideVolVars.velocity(1);

        // Distance between the two adjacent cell centres
        const auto distance = cellCenterToCellCenter.two_norm();

        // Factor that takes the direction of the unit vector into account
        const auto direction = (unitNormal*cellCenterToCellCenter)/distance;

        // Compute the velocity gradients normal to the face
        const auto gradU = (velocityXRight-velocityXLeft)*direction/distance;
        const auto gradV = (velocityYRight-velocityYLeft)*direction/distance;

        // Use a harmonic average of the depth at the interface.
        const auto averageDepth = 2.0*(waterDepthLeft*waterDepthRight)/(waterDepthLeft + waterDepthRight);

        // Now get the (constant) background turbulent viscosity
        static const auto turbBGViscosity = getParam<Scalar>("ShallowWater.TurbulentViscosity", 1.0e-6);

        Scalar turbViscosity = 0.0;

        // constant eddy viscosity equal to the prescribed background eddy viscosity
        if (!useMixingLengthTurbulenceModel)
            turbViscosity = turbBGViscosity;

        // turbulence model based on mixing length
        else
        {
            // Compute the turbulent viscosity using a combined horizonal/vertical mixing length approach
            // Turbulence coefficients: vertical (Elder like) and horizontal (Smagorinsky like)
            static const auto turbConstV = getParam<Scalar>("ShallowWater.VerticalCoefficientOfMixingLengthModel", 1.0);
            static const auto turbConstH = getParam<Scalar>("ShallowWater.HorizontalCoefficientOfMixingLengthModel", 0.1);

            /** The vertical (Elder-like) contribution to the turbulent viscosity scales with water depth \f[ h \f] and shear velocity \f[ u_{*} \f] :
            *
            * \f[
            * \nu_t^v = c^v \frac{\kappa}{6} u_{*} h
            * \f]
            */
            static const Scalar kappa = 0.41;
            const auto turbViscosityV = turbConstV * (kappa/6.0) * ustar * averageDepth;

            /** The horizontal (Smagorinsky-like) contribution to the turbulent viscosity scales with the water depth (squared)
            * and the magnitude of the stress (rate-of-strain) tensor:
            *
            * \f[
            * nu_t^h = (c^h h)^2 \sqrt{ 2\left(\frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)^2 + 2\left(\frac{\partial v}{\partial y}\right)^2 }
            * \f]
            *
            * However, based on the velocity vectors in the direct neighbours of the volume face, it is not possible to compute all components of the stress tensor.
            * Only the gradients of u and v in the direction of the vector between the cell centres is available.
            * To avoid the reconstruction of the full velocity gradient tensor based on a larger stencil,
            * the horizontal contribution to the eddy viscosity (in the mixing-length model) is computed using only the velocity gradients normal to the face:
            *
            * \f[
            * \frac{\partial u}{\partial n} , \frac{\partial v}{\partial n}
            * \f]
            *
            * In other words, the present approximation of the horizontal contribution to the turbulent viscosity reduces to:
            *
            * \f[
            * nu_t^h = (c^h h)^2 \sqrt{ 2\left(\frac{\partial u}{\partial n}\right)^2 + 2\left(\frac{\partial v}{\partial n}\right)^2 }
            * \f]
            *
            It should be noted that this simplified approach is formally inconsistent and will result in a turbulent viscosity that is dependent on the grid (orientation).
            */
            const auto mixingLengthSquared = turbConstH * turbConstH * averageDepth * averageDepth;
            const auto turbViscosityH      = mixingLengthSquared * sqrt(2.0*gradU*gradU + 2.0*gradV*gradV);

            // Total turbulent viscosity
            turbViscosity = turbBGViscosity + sqrt(turbViscosityV*turbViscosityV + turbViscosityH*turbViscosityH);
        }

        // Compute the viscous momentum fluxes
        const auto uViscousFlux = turbViscosity * averageDepth * gradU;
        const auto vViscousFlux = turbViscosity * averageDepth * gradV;

        // compute the mobility of the flux with the fluxlimiter
        static const auto upperWaterDepthFluxLimiting = getParam<double>("FluxLimiterLET.UpperWaterDepth", 1e-3);
        static const auto lowerWaterDepthFluxLimiting = getParam<double>("FluxLimiterLET.LowerWaterDepth", 1e-5);

        const auto limitingDepth = (waterDepthLeft + waterDepthRight) * 0.5;
        const auto mobility = ShallowWater::fluxLimiterLET(limitingDepth,
                                                             limitingDepth,
                                                             upperWaterDepthFluxLimiting,
                                                             lowerWaterDepthFluxLimiting);

        localFlux[0] = 0.0;
        localFlux[1] = -uViscousFlux * mobility * scvf.area();
        localFlux[2] = -vViscousFlux * mobility * scvf.area();

        return localFlux;
    }
};

} // end namespace Dumux

#endif
