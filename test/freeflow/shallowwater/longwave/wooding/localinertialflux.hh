// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Local-inertial advective flux for the shallow water equations
 */
#ifndef DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_LOCALINERTIALFLUX_HH
#define DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_LOCALINERTIALFLUX_HH

#include <cmath>

#include <dumux/common/parameters.hh>
#include <dumux/flux/shallowwaterflux.hh>

namespace Dumux::Wooding {

/*!
 * \ingroup ShallowWaterTests
 * \brief Advective flux that drops momentum advection and drives momentum by the
 *        free-surface gradient
 *
 * This is the local-inertial approximation of Bates, Horritt and Fewtrell (2010): local
 * acceleration, the free-surface gradient and friction are kept, the advective term
 * \f$ \nabla \cdot (\mathbf{q} \otimes \mathbf{u}) \f$ is dropped. It lies between the
 * diffusive wave, which drops local acceleration as well, and the full equations.
 *
 * The hydrostatic reconstruction of the default flux reconstructs a face dry where the bed step
 * between the cells exceeds the depth, as on the planes of this catchment. Here, momentum is
 * driven by \f$ g h \nabla H \f$ in Green-Gauss form, with the face value taken relative to the
 * cell it acts on, which vanishes for water at rest on any bed.
 */
template<class NumEqVector>
class LocalInertialFlux
{
    using DefaultFlux = ShallowWaterFlux<NumEqVector>;

public:
    using Cache = typename DefaultFlux::Cache;
    using CacheFiller = typename DefaultFlux::CacheFiller;

    template<class Problem, class FVElementGeometry, class ElementVolumeVariables>
    static NumEqVector flux(const Problem& problem,
                            const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const auto& nxy = scvf.unitOuterNormal();
        const auto gravity = problem.spatialParams().gravity(scvf.center());

        using std::max, std::min, std::sqrt;
        const auto freeSurfaceInside = insideVolVars.bedSurface() + insideVolVars.waterDepth();
        const auto freeSurfaceOutside = outsideVolVars.bedSurface() + outsideVolVars.waterDepth();

        // depth above the higher bed, with a smooth maximum of the free surfaces, since Newton's
        // method with a numerically differentiated Jacobian does not converge across the kink of
        // a plain maximum
        static const auto eps = getParam<typename NumEqVector::value_type>(
            "ShallowWater.Regularization.FreeSurfaceEpsilon", 1e-4
        );
        const auto surfaceDifference = freeSurfaceOutside - freeSurfaceInside;
        const auto higherFreeSurface = 0.5*(freeSurfaceInside + freeSurfaceOutside)
                                     + 0.5*sqrt(surfaceDifference*surfaceDifference + eps*eps);
        const auto flowDepth = max(0.0, higherFreeSurface
                                        - max(insideVolVars.bedSurface(), outsideVolVars.bedSurface()));

        const auto normalVelocity = 0.5*(
            (insideVolVars.velocity(0) + outsideVolVars.velocity(0))*nxy[0]
          + (insideVolVars.velocity(1) + outsideVolVars.velocity(1))*nxy[1]
        );
        const auto force = gravity*flowDepth*0.5*surfaceDifference;

        const auto upwindDepth = normalVelocity > 0.0 ? insideVolVars.waterDepth()
                                                      : outsideVolVars.waterDepth();

        NumEqVector flux(0.0);
        flux[0] = min(upwindDepth, flowDepth)*normalVelocity*scvf.area();
        flux[1] = force*nxy[0]*scvf.area();
        flux[2] = force*nxy[1]*scvf.area();
        return flux;
    }
};

} // end namespace Dumux::Wooding

#endif
