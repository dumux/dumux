// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \copydoc Dumux::ShallowWaterViscousFlux
 */
#ifndef DUMUX_FLUX_SHALLOW_WATER_VISCOUS_FLUX_HH
#define DUMUX_FLUX_SHALLOW_WATER_VISCOUS_FLUX_HH

#include <cmath>
#include <algorithm>
#include <utility>
#include <type_traits>
#include <array>

#include <dune/common/std/type_traits.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>
#include <dumux/flux/fluxvariablescaching.hh>
#include <dumux/flux/shallowwater/fluxlimiterlet.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {
// helper struct detecting if the user-defined spatial params class has a frictionLaw function
template <typename T, typename ...Ts>
using FrictionLawDetector = decltype(std::declval<T>().frictionLaw(std::declval<Ts>()...));

template<class T, typename ...Args>
static constexpr bool implementsFrictionLaw()
{ return Dune::Std::is_detected<FrictionLawDetector, T, Args...>::value; }
} // end namespace Detail
#endif

/*!
 * \ingroup Flux
 * \brief Compute the shallow water viscous momentum flux due to viscosity.
 *
 * The viscous momentum flux
 * \f[
 * \int \int_{V} \mathbf{\nabla} \cdot \nu_t h \mathbf{\nabla} \mathbf{u} dV
 * \f]
 * is re-written using Gauss' divergence theorem to:
 * \f[
 * \int_{S_f} \nu_t h \mathbf{\nabla} \mathbf{u} \cdot \mathbf{n_f} dS
 * \f]
 *
 * The effective kinematic viscosity \f$ \nu_t \f$
 * can be calculated by adding a vertical (Elder-like)
 * and a horizontal (Smagorinsky-like) part. This enabled by setting
 * "ShallowWater.UseMixingLengthTurbulenceModel" to "true".
 *
 * For now the calculation of the shallow water viscous momentum flux is implemented
 * strictly for 2D depth-averaged models (i.e. 3 equations).
 */
template<class NumEqVector, typename std::enable_if_t<NumEqVector::size() == 3, int> = 0>
class ShallowWaterViscousFlux
{

public:

    using Cache = FluxVariablesCaching::EmptyDiffusionCache;
    using CacheFiller = FluxVariablesCaching::EmptyCacheFiller;
    /*!
     * \ingroup Flux
     * \brief Compute the viscous momentum flux contribution from the interface
     *        shear stress
     */
    template<class Problem, class FVElementGeometry, class ElementVolumeVariables>
    static NumEqVector flux(const Problem& problem,
                            const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {
        using Scalar = typename NumEqVector::value_type;
        using FluidSystem = typename ElementVolumeVariables::VolumeVariables::FluidSystem;

        NumEqVector localFlux(0.0);

        // Get the inside and outside volume variables
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());

        const auto gradVelocity = [&]()
        {
            // The left (inside) and right (outside) states
            const auto velocityXLeft   = insideVolVars.velocity(0);
            const auto velocityYLeft   = insideVolVars.velocity(1);
            const auto velocityXRight  = outsideVolVars.velocity(0);
            const auto velocityYRight  = outsideVolVars.velocity(1);

            // Compute the velocity gradients normal to the face
            // Factor that takes the direction of the unit vector into account
            const auto& cellCenterToCellCenter = outsideScv.center() - insideScv.center();
            const auto distance = cellCenterToCellCenter.two_norm();
            const auto& unitNormal = scvf.unitOuterNormal();
            const auto direction = (unitNormal*cellCenterToCellCenter)/distance;
            return std::array<Scalar, 2>{
                (velocityXRight-velocityXLeft)*direction/distance,
                (velocityYRight-velocityYLeft)*direction/distance
            };
        }();

        // Use a harmonic average of the depth at the interface.
        const auto waterDepthLeft  = insideVolVars.waterDepth();
        const auto waterDepthRight = outsideVolVars.waterDepth();
        const auto averageDepth = 2.0*(waterDepthLeft*waterDepthRight)/(waterDepthLeft + waterDepthRight);

        const Scalar effectiveKinematicViscosity = [&]()
        {
            // The background viscosity, per default this is just the fluid viscosity (e.g. for the viscous flow regime)
            // The default may be overwritten by the user, by setting the parameter "ShallowWater.TurbulentViscosity" to
            // enable setting a background viscosity that already contains some effects of turbulence
            static const auto backgroundKinematicViscosity = getParamFromGroup<Scalar>(
                problem.paramGroup(), "ShallowWater.TurbulentViscosity", insideVolVars.viscosity()/insideVolVars.density()
            );

            // Check whether the mixing-length turbulence model is used
            static const auto useMixingLengthTurbulenceModel = getParamFromGroup<bool>(
                problem.paramGroup(), "ShallowWater.UseMixingLengthTurbulenceModel", false
            );

            // constant eddy viscosity equal to the prescribed background eddy viscosity
            if (!useMixingLengthTurbulenceModel)
                return backgroundKinematicViscosity;

            using SpatialParams = typename Problem::SpatialParams;
            using E = typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity;
            using SCV = typename FVElementGeometry::SubControlVolume;
            if constexpr (!Detail::implementsFrictionLaw<SpatialParams, E, SCV>())
                DUNE_THROW(Dune::IOError, "Mixing length turbulence model enabled but spatial parameters do not implement the frictionLaw interface!");
            else
            {
                // turbulence model based on mixing length
                // Compute the turbulent viscosity using a combined horizontal/vertical mixing length approach
                // Turbulence coefficients: vertical (Elder like) and horizontal (Smagorinsky like)
                static const auto turbConstV = getParamFromGroup<Scalar>(problem.paramGroup(), "ShallowWater.VerticalCoefficientOfMixingLengthModel", 1.0);
                static const auto turbConstH = getParamFromGroup<Scalar>(problem.paramGroup(), "ShallowWater.HorizontalCoefficientOfMixingLengthModel", 0.1);

                /** The vertical (Elder-like) contribution to the turbulent viscosity scales with water depth \f[ h \f] and shear velocity \f[ u_{*} \f] :
                *
                * \f[
                * \nu_t^v = c^v \frac{\kappa}{6} u_{*} h
                * \f]
                */
                constexpr Scalar kappa = 0.41;
                // Compute the average shear velocity on the face
                const Scalar ustar = [&]()
                {
                    // Get the bottom shear stress in the two adjacent cells
                    // Note that element is not needed in spatialParams().frictionLaw (should be removed). For now we simply pass the same element twice
                    const auto& bottomShearStressInside = problem.spatialParams().frictionLaw(element, insideScv).bottomShearStress(insideVolVars);
                    const auto& bottomShearStressOutside = problem.spatialParams().frictionLaw(element, outsideScv).bottomShearStress(outsideVolVars);
                    const auto bottomShearStressInsideMag = bottomShearStressInside.two_norm();
                    const auto bottomShearStressOutsideMag = bottomShearStressOutside.two_norm();

                    // Use a harmonic average of the viscosity and the depth at the interface.
                    using std::max;
                    const auto averageBottomShearStress = 2.0*(bottomShearStressInsideMag*bottomShearStressOutsideMag)
                            / max(1.0e-8,bottomShearStressInsideMag + bottomShearStressOutsideMag);

                    // Shear velocity possibly needed in mixing-length turbulence model in the computation of the diffusion flux
                    using std::sqrt;
                    static_assert(!FluidSystem::isCompressible(0),
                        "The shallow water model assumes incompressible fluids"
                    );
                    // Since the shallow water equations are incompressible, the density is constant.
                    // Therefore, insideVolVars.density() is equal to outsideVolVars.density().
                    return sqrt(averageBottomShearStress/insideVolVars.density());
                }();

                const auto turbViscosityV = turbConstV * (kappa/6.0) * ustar * averageDepth;

                /** The horizontal (Smagorinsky-like) contribution to the turbulent viscosity scales with the water depth (squared)
                * and the magnitude of the stress (rate-of-strain) tensor:
                *
                * \f[
                * \nu_t^h = (c^h h)^2 \sqrt{ 2\left(\frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)^2 + 2\left(\frac{\partial v}{\partial y}\right)^2 }
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
                * \nu_t^h = (c^h h)^2 \sqrt{ 2\left(\frac{\partial u}{\partial n}\right)^2 + 2\left(\frac{\partial v}{\partial n}\right)^2 }
                * \f]
                *
                It should be noted that this simplified approach is formally inconsistent and will result in a turbulent viscosity that is dependent on the grid (orientation).
                */
                using std::sqrt;
                const auto [gradU, gradV] = gradVelocity;
                const auto mixingLengthSquared = turbConstH * turbConstH * averageDepth * averageDepth;
                const auto turbViscosityH = mixingLengthSquared * sqrt(2.0*gradU*gradU + 2.0*gradV*gradV);

                // Total turbulent viscosity
                return backgroundKinematicViscosity + sqrt(turbViscosityV*turbViscosityV + turbViscosityH*turbViscosityH);
            }
        }();

        // Compute the viscous momentum fluxes with connecting water depth at the interface
        using std::min;
        using std::max;
        const auto freeSurfaceInside = insideVolVars.waterDepth() + insideVolVars.bedSurface();
        const auto freeSurfaceOutside = outsideVolVars.waterDepth() + outsideVolVars.bedSurface();
        const auto interfaceWaterDepth = max(min(freeSurfaceInside , freeSurfaceOutside) - max(insideVolVars.bedSurface(),outsideVolVars.bedSurface()),0.0);
        const auto& [gradU, gradV] = gradVelocity;
        const auto uViscousFlux = effectiveKinematicViscosity * interfaceWaterDepth * gradU;
        const auto vViscousFlux = effectiveKinematicViscosity * interfaceWaterDepth * gradV;

        localFlux[0] = 0.0;
        localFlux[1] = -uViscousFlux * scvf.area();
        localFlux[2] = -vViscousFlux * scvf.area();

        return localFlux;
    }
};

} // end namespace Dumux

#endif
